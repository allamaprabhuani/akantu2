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
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "distributed_synchronizer.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_triangle_6;
#endif //AKANTU_USE_IOHELPER
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::initialize(argc, argv);

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator * comm = 
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  // std::stringstream filename;
  // filename << "log-" << prank << ".txt";
  // akantu::debug::setLogFile(filename.str());
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
#ifdef AKANTU_USE_IOHELPER
  akantu::UInt nb_element = mesh.getNbElement(type);
#endif //AKANTU_USE_IOHELPER

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

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(mesh.getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d_para");
  dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
			 iohelper::TRIANGLE2, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getMass().values,
			  spatial_dimension, "mass");
  dumper.AddNodeDataField(model.getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model.getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  akantu::UInt  nb_quadrature_points = model.getFEM().getNbQuadraturePoints(type);
  double * part = new double[nb_element*nb_quadrature_points];
  for (unsigned int i = 0; i < nb_element; ++i)
    for (unsigned int q = 0; q < nb_quadrature_points; ++q)
      part[i*nb_quadrature_points + q] = prank;

  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

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

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) {
      dumper.Dump();
    }
#endif //AKANTU_USE_IOHELPER
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

#ifdef AKANTU_USE_IOHELPER
  delete [] part;
#endif //AKANTU_USE_IOHELPER

  akantu::finalize();
  return EXIT_SUCCESS;
}
