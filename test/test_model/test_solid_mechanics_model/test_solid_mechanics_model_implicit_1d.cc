/**
 * @file   test_solid_mechanics_model_implicit_1d.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Feb 21 14:56:16 2011
 *
 * @brief  test of traction in implicit
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
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

int main(int argc, char *argv[])
{

  akantu::initialize(&argc, &argv);

#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_segment_2;
  akantu::UInt paraview_type = LINE1;
#endif //AKANTU_USE_IOHELPER
  akantu::UInt spatial_dimension = 1;
  // akantu::UInt max_steps = 10000;
  // akantu::Real time_factor = 0.2;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("segment1.msh", mesh);

// #ifdef AKANTU_USE_SCOTCH
//   mesh_io.write("bar1_breorder.msh", mesh);

//   akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
//   partition->reorder();
//   delete partition;

//   mesh_io.write("bar1_reorder.msh", mesh);
// #endif

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);
#endif //AKANTU_USE_IOHELPER

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

  model->readMaterials("material.dat");
  model->initMaterials();

  model->initImplicit();

  std::cout << model->getMaterial(0) << std::endl;

  /// boundary conditions
  // akantu::Real eps = 1e-16;
  // for (akantu::UInt i = 0; i < nb_nodes; ++i) {
  //   if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= (bar_length  - bar_length / 10.))
  //     model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - (bar_length - bar_length / 10.)) / 100.;

  //   if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
  // 	model->getBoundary().values[spatial_dimension*i] = true;

  //   if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
  //      model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
  //     model->getBoundary().values[spatial_dimension*i + 1] = true;
  //   }
  // }

  //  model->getForce().at(0) = -1000;
  model->getBoundary().at(0,0) = true;
  model->getForce().at(1,0) = 1000;

  model->initializeUpdateResidualData();

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "implicit");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "forces");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetEmbeddedValue("forces", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  //  dumper.Dump();
#endif //AKANTU_USE_IOHELPER


  akantu::debug::setDebugLevel(akantu::dblInfo);
  akantu::UInt count = 0;
  model->updateResidual();
#ifdef AKANTU_USE_IOHELPER
    dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  while(!model->testConvergenceResidual(1e-1) && (count <= 10)) {
    std::cout << "Iter : " << ++count << std::endl;
    model->assembleStiffnessMatrix();

    model->solveStatic();

    model->getStiffnessMatrix().saveMatrix("Ktmp.mtx");

    model->updateResidual();
    //    model->assembleStiffnessMatrix();
#ifdef AKANTU_USE_IOHELPER
    dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  delete model;
  akantu::finalize();

  return EXIT_SUCCESS;
}
