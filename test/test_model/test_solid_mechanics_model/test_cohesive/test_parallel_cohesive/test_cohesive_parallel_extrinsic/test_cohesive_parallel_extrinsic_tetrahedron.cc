/**
 * @file   test_cohesive_parallel_extrinsic_tetrahedron.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Sep 20 16:53:10 2013
 *
 * @brief  3D extrinsic cohesive elements test
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "dumper_paraview.hh"
#include "static_communicator.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);
  debug::setDebugLevel(dblWarning);

  const UInt max_steps = 1000;

  UInt spatial_dimension = 3;
  Mesh mesh(spatial_dimension);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("tetrahedron.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    //    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    //    debug::setDebugLevel(dblWarning);
  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition, NULL, true);

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh << std::endl;
  // debug::setDebugLevel(dblWarning);


  model.initFull("material.dat", _explicit_lumped_mass, true);

  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  Array<Real> limits(spatial_dimension, 2);
  limits(0, 0) = -0.01;
  limits(0, 1) = 0.01;
  limits(1, 0) = -100;
  limits(1, 1) = 100;
  limits(2, 0) = -100;
  limits(2, 1) = 100;

  model.enableFacetsCheckOnArea(limits);

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  // debug::setDebugLevel(dblDump);
  // std::cout << mesh_facets << std::endl;
  // debug::setDebugLevel(dblWarning);

  Real time_step = model.getStableTimeStep()*0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  /// initial conditions
  Real loading_rate = 0.5;
  Real disp_update = loading_rate * time_step;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  model.setBaseName("extrinsic_parallel_tetrahedron");
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("velocity"    );
  model.addDumpFieldVector("acceleration");
  model.addDumpFieldVector("residual"    );
  model.addDumpFieldTensor("stress");
  model.addDumpFieldTensor("strain");
  model.addDumpField("partitions");
  //  model.getDumper().getDumper().setMode(iohelper::BASE64);
  model.dump();

  DumperParaview dumper("cohesive_elements");
  dumper.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_cohesive);
  DumperIOHelper::Field * cohesive_displacement =
    new DumperIOHelper::NodalField<Real>(model.getDisplacement());
  cohesive_displacement->setPadding(3);
  dumper.registerField("displacement", cohesive_displacement);
  // dumper.registerField("damage", new DumperIOHelper::
  // 		       HomogenizedField<Real,
  // 					DumperIOHelper::InternalMaterialField>(model,
  // 									       "damage",
  // 									       spatial_dimension,
  // 									       _not_ghost,
  // 									       _ek_cohesive));
  dumper.dump();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
	displacement(n, 0) += disp_update * position(n, 0);
    }

    model.checkCohesiveStress();

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 10 == 0) {
      model.dump();
      dumper.dump();
      if(prank == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 0) + displacement(n, 0) > 0) {
    // 	displacement(n, 0) -= 0.01;
    //   }
    // }

    //    Real Ed = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getDissipatedEnergy();
    //    Real Er = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;

  }

  dumper.dump();

  // edis.close();
  // erev.close();

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 400;

  if(prank == 0) {
    std::cout << Ed << " " << Edt << std::endl;

    if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
      std::cout << "The dissipated energy is incorrect" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
