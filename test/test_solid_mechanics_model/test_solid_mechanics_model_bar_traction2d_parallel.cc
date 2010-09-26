/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * <insert license here>
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
#include "communicator.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

#define CHECK_STRESS

int main(int argc, char *argv[])
{
  akantu::ElementType type = akantu::_triangle_1;
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 1000;
  akantu::Real time_factor = 0.2;

  akantu::initialize(&argc, &argv);

  akantu::debug::setDebugLevel(akantu::dblWarning);
  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::Communicator * communicator;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("bar.msh", mesh);
    akantu::MeshPartition * partition =
      new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  //  std::cout << mesh << std::endl;

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);

  /// model initialization
  model->initVectors();

  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));


  model->readMaterials("material.dat");
  model->initMaterials();
  model->registerSynchronizer(*communicator);

  model->initModel();

  if(prank == 0)
    std::cout << model->getMaterial(0) << std::endl;

  model->assembleMass();


  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= 9)
      model->getDisplacement().values[spatial_dimension*i]
	= (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - 9) / 100.;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
      model->getBoundary().values[spatial_dimension*i] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  model->synchronizeBoundaries();


  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d_para");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 TRIANGLE1, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getMass().values,
			  1, "mass");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  //dumper.Dump();

  DumperParaview dumper_ghost;
  if (psize > 1) {
    unsigned int nb_ghost_element = mesh.getNbGhostElement(type);
    dumper_ghost.SetMode(TEXT);
    dumper_ghost.SetParallelContext(prank, psize);
    dumper_ghost.SetPoints(mesh.getNodes().values,
			   spatial_dimension, nb_nodes, "bar2d_para_ghost");
    dumper_ghost.SetConnectivity((int*) mesh.getGhostConnectivity(type).values,
				 TRIANGLE1, nb_ghost_element, C_MODE);
    dumper_ghost.AddNodeDataField(model->getDisplacement().values,
				  spatial_dimension, "displacements");
    dumper_ghost.AddNodeDataField(model->getMass().values,
				  1, "mass");
    dumper_ghost.AddNodeDataField(model->getVelocity().values,
				  spatial_dimension, "velocity");
    dumper_ghost.AddNodeDataField(model->getResidual().values,
				  spatial_dimension, "force");
    dumper_ghost.AddNodeDataField(model->getAcceleration().values,
				  spatial_dimension, "acceleration");
    dumper_ghost.SetPrefix("paraview/");
    dumper_ghost.Init();
    //    dumper_ghost.Dump();
  }
#endif //AKANTU_USE_IOHELPER

  model->setPotentialEnergyFlagOn();

  std::ofstream energy;
  if(prank == 0) {
    energy.open("energy.csv");
    energy << "id,epot,ekin,tot" << std::endl;
  }

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    if(prank == 0 && (s % 10 == 0)) {
      energy << s << "," << epot << "," << ekin << "," << epot + ekin
	     << std::endl;
    }

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) {
      dumper.Dump();
      //      if (psize > 1) dumper_ghost.Dump();
    }
#endif //AKANTU_USE_IOHELPER
    if(prank == 0)
      if(s % 10 == 0) std::cerr << "passing step " << s << "/" << max_steps << std::endl;
  }

  delete [] part;
  delete model;

  akantu::finalize();
  return EXIT_SUCCESS;
}
