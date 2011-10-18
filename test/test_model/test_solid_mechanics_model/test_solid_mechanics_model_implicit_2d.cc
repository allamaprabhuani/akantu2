/**
 * @file   test_solid_mechanics_model_implicit_2d.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Feb 14 14:56:16 2011
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
#include "static_communicator.hh"
#include "distributed_synchronizer.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

#define bar_length 1
#define bar_height 1

// static void traction(__attribute__ ((unused)) double * position,double * stress){
//    memset(stress,0,sizeof(akantu::Real)*4);
//    if((fabs(position[1] - bar_height) < akantu::Math::tolerance) || (fabs(position[0] - bar_length) < akantu::Math::tolerance)) {
//      stress[0] = 1000;
//      stress[3] = 1000;
//    }
// }

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(&argc, &argv);

  akantu::UInt spatial_dimension = 2;

  akantu::Mesh mesh(spatial_dimension);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("square_implicit2.msh", mesh);

    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    //   partition->reorder();
    partition->partitionate(psize);
  }
  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);
  model->initParallel(partition);
  delete partition;

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
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

  if (prank == 0)
    std::cout << model->getMaterial(0) << std::endl;

  /// boundary conditions
  const  akantu::Vector<akantu::Real> & position = mesh.getNodes();
  akantu::Vector<bool> & boundary = model->getBoundary();
  akantu::Vector<akantu::Real> & displacment = model->getDisplacement();

  for (akantu::UInt n = 0; n < nb_nodes; ++n) {
    if(position(n,0) < akantu::Math::getTolerance()) boundary(n,0) = true;
    if(position(n,1) < akantu::Math::getTolerance()) boundary(n,1) = true;

    if(std::abs(position(n,0) - bar_length) < akantu::Math::getTolerance()) {
      boundary(n,0) = true;
      displacment(n,0) = 0.1;
    }
  }

#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_triangle_6;
  akantu::UInt paraview_type = TRIANGLE2;
  akantu::UInt nb_element = model->getFEM().getMesh().getNbElement(type);

  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetParallelContext(prank, psize);
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
#endif //AKANTU_USE_IOHELPER


  akantu::UInt count = 0;
  model->updateResidual();

#ifdef AKANTU_USE_IOHELPER
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  akantu::Real norm;
  model->assembleStiffnessMatrix();
  while(!model->testConvergenceResidual(1e-3, norm) && (count < 100)) {
    if (prank == 0)
      std::cout << "Iter : " << ++count << " - residual norm : " << norm << std::endl;

    model->solveStatic();
    model->updateResidual();
  };

#ifdef AKANTU_USE_IOHELPER
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  delete model;
  akantu::finalize();

  return EXIT_SUCCESS;
}
