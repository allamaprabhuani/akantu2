/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "mesh_partition_scotch.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_iohelper.hh"
#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::initialize(argc, argv);

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator & comm =
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm.getNbProc();
  akantu::Int prank = comm.whoAmI();

  akantu::debug::setDebugLevel(akantu::dblWarning);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("bar2.msh", mesh);
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  akantu::SolidMechanicsModel model(mesh);

  /// model initialization
  model.initParallel(partition);

  model.initFull("material.dat");

  if(prank == 0)
    std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();


  akantu::UInt nb_nodes = mesh.getNbNodes();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  akantu::Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & pos = mesh.getNodes();
  akantu::Vector<akantu::Real> & disp = model.getDisplacement();
  akantu::Vector<bool> & boun = model.getBoundary();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  model.synchronizeBoundaries();

  model.updateResidual();

  model.setBaseName("bar2d_parallel");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.addDumpField("partitions"  );
  model.dump();

  std::ofstream energy;
  if(prank == 0) {
    energy.open("energy_bar_2d_para.csv");
    energy << "id,rtime,epot,ekin,tot" << std::endl;
  }

  double total_time = 0.;

  /// Setting time step
  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    double start = MPI_Wtime();

    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    double end = MPI_Wtime();

    akantu::Real epot = model.getPotentialEnergy();
    akantu::Real ekin = model.getKineticEnergy();

    total_time += end - start;

    if(prank == 0 && (s % 100 == 0)) {
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    }
    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	     << std::endl;

    if(s % 100 == 0) {
      model.dump();
    }
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

  akantu::finalize();
  return EXIT_SUCCESS;
}
