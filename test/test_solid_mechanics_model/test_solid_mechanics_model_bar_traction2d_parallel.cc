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

int main(int argc, char *argv[])
{
  akantu::ElementType type = akantu::_triangle_2;
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::initialize(&argc, &argv);

  akantu::debug::setDebugLevel(akantu::dblWarning);
  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  akantu::Communicator * communicator;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("bar2.msh", mesh);
    akantu::MeshPartition * partition =
      new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator = akantu::Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  /// model initialization
  model->initVectors();

  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getResidual().values,        0,
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

  model->assembleMassLumped();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    // model->getDisplacement().values[spatial_dimension*i]
    //   = model->getFEM().getMesh().getNodes().values[spatial_dimension*i] / 100.;
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


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d_para");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 TRIANGLE2, nb_element, C_MODE);
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
  akantu::UInt  nb_quadrature_points = akantu::FEM::getNbQuadraturePoints(type);
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

  model->setPotentialEnergyFlagOn();
  std::ofstream energy;
  if(prank == 0) {
    energy.open("energy.csv");
    energy << "id,epot,ekin,tot" << std::endl;
  }

  double total_time = 0.;

  /// Setting time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  if(prank == 0)
    std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    double start = MPI_Wtime();

    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    double end = MPI_Wtime();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    total_time += end - start;

    if(prank == 0 && (s % 10 == 0)) {
       std::cerr << "passing step " << s << "/" << max_steps << std::endl;
       energy << s << "," << epot << "," << ekin << "," << epot + ekin
	      << std::endl;
    }

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) {
      dumper.Dump();
    }
#endif //AKANTU_USE_IOHELPER
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

  delete [] part;
  delete model;

  akantu::finalize();
  return EXIT_SUCCESS;
}
