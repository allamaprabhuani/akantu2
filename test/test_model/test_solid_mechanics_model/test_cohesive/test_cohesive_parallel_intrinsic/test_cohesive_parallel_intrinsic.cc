/**
 * @file   test_cohesive_parallel_intrinsic.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Wed Nov 28 16:59:11 2012
 *
 * @brief  parallel test for intrinsic cohesive elements
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
#include "model.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "dumper_paraview.hh"
#include "static_communicator.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);
  debug::setDebugLevel(dblInfo);

  const UInt max_steps = 350;

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");

  ElementType type = _triangle_6;
  ElementType type_facet = Mesh::getFacetType(type);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    /// insert cohesive elements
    MeshUtils::buildAllFacets(mesh, mesh_facets);

    UInt nb_facet = mesh_facets.getNbElement(type_facet);

    Vector<UInt> facet_insertion;
    Real * bary_facet = new Real[spatial_dimension];
    for (UInt f = 0; f < nb_facet; ++f) {
      mesh_facets.getBarycenter(f, type_facet, bary_facet);
      if (bary_facet[0] > -0.26 && bary_facet[0] < -0.24) facet_insertion.push_back(f);
    }
    delete[] bary_facet;

    MeshUtils::insertIntrinsicCohesiveElements(mesh,
    					       mesh_facets,
    					       type_facet,
    					       facet_insertion);


    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    debug::setDebugLevel(dblDump);
    partition->partitionate(psize);
    debug::setDebugLevel(dblInfo);
  }

  SolidMechanicsModelCohesive model(mesh);

  model.initParallel(partition);

  debug::setDebugLevel(dblDump);
  std::cout << mesh << std::endl;
  debug::setDebugLevel(dblInfo);


  model.initFull("material.dat");

  Real time_step = model.getStableTimeStep()*0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  Vector<Real> & position = mesh.getNodes();
  Vector<Real> & velocity = model.getVelocity();
  Vector<bool> & boundary = model.getBoundary();
  //  Vector<Real> & displacement = model.getDisplacement();
  //  const Vector<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();
  Real epsilon = std::numeric_limits<Real>::epsilon();

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (std::abs(position(n, 0) - 1.) < epsilon)
      boundary(n, 0) = true;
  }

  model.synchronizeBoundaries();
  model.updateResidual();

  model.setBaseName("intrinsic_parallel");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  //model.addDumpField("partitions");
  model.addDumpField("force");
  model.getDumper().getDumper().setMode(iohelper::BASE64);
  model.dump();

  /// initial conditions
  Real loading_rate = .2;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 1 == 0) {
      model.dump();
      if(prank == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) + displacement(n, 1) > 0) {
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

  // edis.close();
  // erev.close();

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 2 * sqrt(2);


  if(prank == 0) {
    std::cout << Ed << " " << Edt << std::endl;

    if (std::abs((Ed - Edt) / Edt) > 0.01 || std::isnan(Ed)) {
      std::cout << "The dissipated energy is incorrect" << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();
  if(prank == 0) std::cout << "OK: Test passed!" << std::endl;
  return EXIT_SUCCESS;
}
