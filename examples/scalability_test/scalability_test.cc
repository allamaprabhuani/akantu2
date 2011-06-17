/**
 * @file   scalability_test.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Feb 22 09:35:58 2011
 *
 * @brief  Test de scalability
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
// #ifdef AKANTU_USE_IOHELPER
// #  include "io_helper.h"
// #endif //AKANTU_USE_IOHELPER


using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  initialize(&argc, &argv);
  StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm->getNbProc();
  Int prank = comm->whoAmI();

  /* -------------------------------------------------------------------------- */

  UInt spatial_dimension = 2;
  ElementType type = _quadrangle_4;
  UInt max_steps = 100;
  Real time_factor = 0.2;

  UInt nex = 100, ney = 100 * psize;
  Real height = 1., width = 1. * psize;
  if(argc == 3 || argc == 4) {
    nex = atoi(argv[1]);
    ney = atoi(argv[2]);
    width = ney / 100.;
    if(argc == 4) {
      max_steps = atoi(argv[3]);
    }
  } else if (argc != 1) {
    std::cout << "Usage : " << argv[0] << " [nb_x (default 100) nb_y (default 100 * nb proc)] [nb_step (default 100)]" << std::endl;
    exit(EXIT_FAILURE);
  }


  /* ------------------------------------------------------------------------ */

  // Real epot, ekin;

  Mesh mesh(spatial_dimension);

  if(prank == 0) {
    std::cout << "Generating mesh..." << std::endl;
    Real height_el = height / nex;
    Real width_el = width / ney;
    UInt nnx = nex + 1, nny = ney + 1;

    Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());
    nodes.resize(nnx * nny);

    mesh.addConnecticityType(type);
    Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(mesh.getConnectivity(type));
    connectivity.resize(nex * ney);

    for (UInt i = 0; i < nny; ++i) {
      for (UInt j = 0; j < nnx; ++j) {
	UInt n = i * nnx + j;
	nodes.at(n, 0) = i * width_el;
	nodes.at(n, 1) = j * height_el;
      }
    }

    for (UInt i = 0; i < ney; ++i) {
      for (UInt j = 0; j < nex; ++j) {
	UInt e = i * nex + j;
	connectivity.at(e, 0) = i * nnx + j;
	connectivity.at(e, 1) = (i + 1) * nnx + j;
	connectivity.at(e, 2) = (i + 1) * nnx + (j + 1);
	connectivity.at(e, 3) = i * nnx + j + 1;
      }
    }

    akantu::MeshIOMSH mesh_io;
    mesh_io.write("bar.msh", mesh);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    std::cout << "Partitioning mesh..." << std::endl;
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }
  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);
  model->initParallel(partition);

  /// model initialization
  model->initVectors();

  /// set vectors to 0
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));


  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();
  model->assembleMassLumped();

  //  std::cout << model->getMaterial(0) << std::endl;

  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= (width * .9))
      model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - .9 * width) / 100. ;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
	model->getBoundary().values[spatial_dimension*i] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= (height - eps) ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

// #ifdef AKANTU_USE_IOHELPER
//   akantu::UInt paraview_type = QUAD1;

//   double * part;
//   akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);

//   /// set to 0 only for the first paraview dump
//   memset(model->getResidual().values, 0,
// 	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
//   memset(model->getMaterial(0).getStrain(type).values, 0,
// 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
//   memset(model->getMaterial(0).getStress(type).values, 0,
// 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));

//   DumperParaview dumper;
//   dumper.SetMode(BASE64);
//   dumper.SetParallelContext(prank, psize);
//   dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
// 		   spatial_dimension, nb_nodes, "coordinates");
//   dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
// 			 paraview_type, nb_element, C_MODE);
//   dumper.AddNodeDataField(model->getDisplacement().values,
// 			  spatial_dimension, "displacements");
//   dumper.AddNodeDataField(model->getVelocity().values,
// 			  spatial_dimension, "velocity");
//   dumper.AddNodeDataField(model->getResidual().values,
// 			  spatial_dimension, "force");
//   dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
// 			  spatial_dimension*spatial_dimension, "strain");
//   dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
// 			  spatial_dimension*spatial_dimension, "stress");

//   akantu::UInt  nb_quadrature_points = model->getFEM().getNbQuadraturePoints(type);
//   part = new double[nb_element*nb_quadrature_points];
//   for (unsigned int i = 0; i < nb_element; ++i)
//     for (unsigned int q = 0; q < nb_quadrature_points; ++q)
//       part[i*nb_quadrature_points + q] = prank;

//   dumper.AddElemDataField(part, 1, "partitions");
//   dumper.SetEmbeddedValue("displacements", 1);
//   dumper.SetPrefix("paraview/");
//   dumper.Init();
//   dumper.Dump();
// #endif //AKANTU_USE_IOHELPER

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  model->setTimeStep(time_step);

  // std::ofstream energy;
  // energy.open("energy.csv");
  // energy << "id,epot,ekin,tot" << std::endl;
  double total_time = 0.;
  if(prank == 0) {
    std::cout << "Time Step = " << time_step << "s" << std::endl;
    std::cerr << "passing step " << 0 << "/" << max_steps << std::endl;
  }
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    double start = MPI_Wtime();

    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    double end = MPI_Wtime();

    total_time += end - start;
    // epot = model->getPotentialEnergy();
    // ekin = model->getKineticEnergy();

    if(s % (std::max(1,(int)max_steps/10)) == 0 && prank == 0) {
      std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    }
    // energy << s << "," << epot << "," << ekin << "," << epot + ekin
    // 	   << std::endl;

// #ifdef AKANTU_USE_IOHELPER
//     if(s % 10 == 0) {
//       dumper.Dump();
//     }
// #endif //AKANTU_USE_IOHELPER
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

// #ifdef AKANTU_USE_IOHELPER
//   delete [] part;
// #endif //AKANTU_USE_IOHELPER
  // energy.close();

  finalize();


  return EXIT_SUCCESS;
}
