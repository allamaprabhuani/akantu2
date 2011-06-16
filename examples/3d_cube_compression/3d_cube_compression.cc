/**
 * @file   3d_cube_compression.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 24 23:15:47 2010
 *
 * @brief  3d cube compression
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
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER


int main(int argc, char *argv[])
{
  akantu::ElementType type = akantu::_tetrahedron_4;
#ifdef AKANTU_USE_IOHELPER
  akantu::UInt paraview_type = TETRA1;
#endif //AKANTU_USE_IOHELPER
  akantu::UInt spatial_dimension = 3;
  akantu::UInt max_steps = 100;
  akantu::Real time_step = 1e-6;

  akantu::initialize(&argc, &argv);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  //  akantu::Real epot, ekin;
  akantu::Mesh mesh(spatial_dimension);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("cube.msh", mesh);
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }
  

  // akantu::Mesh mesh(spatial_dimension);
  // akantu::MeshIOMSH mesh_io;
  // mesh_io.read("cube.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);
  model->initParallel(partition);

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);

  std::cout << "Nb nodes : " << nb_nodes << " - nb elements : " << nb_element << std::endl;

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
 model->initModel();

  model->initExplicit();
  model->readMaterials("material.dat");
  model->initMaterials();


  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();


#ifdef AKANTU_USE_IOHELPER
  /// set to 0 only for the first paraview dump
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getMaterial(0).getStrain(type).values, 0,
	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
  memset(model->getMaterial(0).getStress(type).values, 0,
	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
#endif //AKANTU_USE_IOHELPER


  /// boundary conditions
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    // if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= 9)
    //   model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - 9) / 100. ;
    if(fabs(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 2] - 1) <= eps)
	model->getForce().values[spatial_dimension*i + 2] = -250;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 2] <= eps)
	model->getBoundary().values[spatial_dimension*i + 2] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 0] <= eps) {
      model->getBoundary().values[spatial_dimension*i + 0] = true;
    }
  }


  //  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);
  //  model->setTimeStep(3.54379e-07);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  double total_time = 0.;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    double start = MPI_Wtime();

    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    double end = MPI_Wtime();
    total_time += end - start;

#ifdef AKANTU_USE_IOHELPER
    if(s % 1 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  if(prank == 0) std::cout << "Time : " << psize << " " << total_time / max_steps << " " << total_time << std::endl;

  akantu::finalize();

  return EXIT_SUCCESS;
}
