/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
  akantu::UInt spatial_dimension = 1;
  akantu::ElementType type = akantu::_line_2;
  akantu::UInt max_steps = 10000;
  akantu::Real time_factor = 0.2;

  akantu::initialize(&argc, &argv);

  //  akantu::debug::setDebugLevel(akantu::dblDump);
  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::Communicator * communicator;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("line.msh", mesh);
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

  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();


#ifdef AKANTU_USE_IOHELPER
  /// set to 0 only for the first paraview dump
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
#endif //AKANTU_USE_IOHELPER


  /// boundary conditions
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    model->getDisplacement().values[spatial_dimension*i] = model->getFEM().getMesh().getNodes().values[i] / 100. ;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= 1e-15)
      model->getBoundary().values[i] = true;
  }


  memset(model->getForce().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));

  //  model->synchronizeBoundaries();


  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "line_para");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 LINE2, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model->getMass().values,
			  spatial_dimension, "mass");
  double * part = new double[nb_element];
  for (unsigned int i = 0; i < nb_element; ++i)
    part[i] = prank;
  dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER


  model->setPotentialEnergyFlagOn();
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    if(prank == 0) {
      std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
		<< std::endl;
    }

#ifdef AKANTU_USE_IOHELPER
    //    if(s % 100 == 0)
    dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cerr << "passing step " << s << "/" << max_steps << std::endl;
  }

  delete [] part;

  akantu::finalize();

  return EXIT_SUCCESS;
}
