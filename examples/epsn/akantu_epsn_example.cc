/**
 * @file   akantu_epsn_example.cc
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
#include <EpsnSimulation_Interface.hh>
#include <RedSYM.hh>

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
#undef AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
  std::string simulation_name = "bar_traction";
  std::string xml_filename = "bar_traction.xml";

  akantu::ElementType type = akantu::_quadrangle_4;
  RedSYM::CellType redsym_type = RedSYM::_quadrilateral;

    //  akantu::ElementType type = akantu::_triangle_3;
    //  RedSYM::CellType redsym_type = RedSYM::_triangle

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
    mesh_io.read("bar.msh", mesh);
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


  EPSN::Simulation::Interface * epsn_itfc;
  epsn_itfc = new EPSN::Simulation::Interface;
  epsn_itfc->initORB(argc,argv);

  if(prank == 0)
    epsn_itfc->initProxyAndNode(simulation_name, xml_filename, prank, psize);
  else
    epsn_itfc->initNode(simulation_name, prank, psize);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  /// model initialization
  model->initVectors();

  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getResidual().values,     0,
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

  //  akantu::UInt  nb_quadrature_points = akantu::FEM::getNbQuadraturePoints(type);

  // double * part = new double[nb_element*nb_quadrature_points];
  // for (unsigned int i = 0; i < nb_element; ++i)
  //   for (unsigned int q = 0; q < nb_quadrature_points; ++q)
  //     part[i*nb_quadrature_points + q] = prank;

  // dumper.AddElemDataField(part, 1, "partitions");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  akantu::UInt nb_nodes_per_element = akantu::Mesh::getNbNodesPerElement(type);
  RedSYM::Mesh * rs_mesh = new RedSYM::Mesh("bar",
					    prank, psize,
					    spatial_dimension,
					    redsym_type,
					    RedSYM::_c_numbering);

  int cell_variable   = rs_mesh->addVariable("connectivity", RedSYM::_long,  nb_nodes_per_element);
  int node_variable   = rs_mesh->addVariable("nodes", RedSYM::_double, spatial_dimension);
  int displacements_variable = rs_mesh->addVariable("displacements", RedSYM::_double, spatial_dimension                  );
  int mass_variable          = rs_mesh->addVariable("mass"         , RedSYM::_double, 1                                  );
  int velocity_variable      = rs_mesh->addVariable("velocity"     , RedSYM::_double, spatial_dimension                  );
  int force_variable         = rs_mesh->addVariable("force"        , RedSYM::_double, spatial_dimension                  );
  int acceleration_variable  = rs_mesh->addVariable("acceleration" , RedSYM::_double, spatial_dimension                  );
  int strain_variable        = rs_mesh->addVariable("strain"       , RedSYM::_double, spatial_dimension*spatial_dimension);
  int stress_variable        = rs_mesh->addVariable("stress"       , RedSYM::_double, spatial_dimension*spatial_dimension);

  rs_mesh->setCellVariable(cell_variable);
  rs_mesh->setNodeVariable(node_variable);

  rs_mesh->setNodeDataVariable(displacements_variable);
  rs_mesh->setNodeDataVariable(mass_variable         );
  rs_mesh->setNodeDataVariable(velocity_variable     );
  rs_mesh->setNodeDataVariable(force_variable        );
  rs_mesh->setNodeDataVariable(acceleration_variable );

  rs_mesh->setCellDataVariable(strain_variable       );
  rs_mesh->setCellDataVariable(stress_variable       );

  RedSYM::Key region = rs_mesh->addRegion(nb_element, nb_nodes);
  rs_mesh->setRegionTag(region, prank*100);
  rs_mesh->setNumberOfCells(region, nb_element);
  rs_mesh->setNumberOfNodes(region, nb_nodes);
  rs_mesh->wrapCells(region, model->getFEM().getMesh().getConnectivity(type).values);
  rs_mesh->wrapNodes(region, model->getFEM().getMesh().getNodes().values);

  rs_mesh->wrap(displacements_variable, region, model->getDisplacement().values             );
  rs_mesh->wrap(mass_variable         , region, model->getMass().values                     );
  rs_mesh->wrap(velocity_variable     , region, model->getVelocity().values                 );
  rs_mesh->wrap(force_variable        , region, model->getResidual().values                 );
  rs_mesh->wrap(acceleration_variable , region, model->getAcceleration().values             );
  rs_mesh->wrap(strain_variable       , region, model->getMaterial(0).getStrain(type).values);
  rs_mesh->wrap(stress_variable       , region, model->getMaterial(0).getStress(type).values);

  epsn_itfc->addData(rs_mesh);

  epsn_itfc->ready();



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
  epsn_itfc->beginHTM();
  epsn_itfc->beginLoop("main");
  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    epsn_itfc->beginTask("body");
    double start = MPI_Wtime();

    epsn_itfc->beginTask("predictor");
    model->explicitPred();
    epsn_itfc->endTask("predictor");

    epsn_itfc->beginTask("updateResidual");
    model->updateResidual();
    epsn_itfc->endTask("updateResidual");

    epsn_itfc->beginTask("updateAcceleration");
    model->updateAcceleration();
    epsn_itfc->endTask("updateAcceleration");

    epsn_itfc->beginTask("corrector");
    model->explicitCorr();
    epsn_itfc->endTask("corrector");

    epsn_itfc->endTask("body");
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

  epsn_itfc->endLoop("main");
  epsn_itfc->endHTM();

  epsn_itfc->finalize();
  epsn_itfc->killORB();


  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

  //  delete [] part;
  delete model;
  energy.close();
  akantu::finalize();
  return EXIT_SUCCESS;
}
