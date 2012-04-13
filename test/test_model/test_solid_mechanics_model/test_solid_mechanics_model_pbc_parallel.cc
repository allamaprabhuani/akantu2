/**
 * @file   test_solid_mechanics_model_pbc_parallel.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Apr 13 16:13:14 2012
 *
 * @brief  test if pbc works in parallel if partition is strips
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_IOHELPER
akantu::ElementType type = akantu::_quadrangle_4;
iohelper::ElemType paraview_type = iohelper::QUAD1;
#endif //AKANTU_USE_IOHELPER

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(iohelper::Dumper & dumper, const akantu::SolidMechanicsModel & model);
static void paraviewDump(iohelper::Dumper & dumper);
#endif

akantu::Vector<akantu::Real> proc_rank(0,1);

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(argc, argv);

  akantu::StaticCommunicator * comm = 
    akantu::StaticCommunicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 1;
  akantu::Real time_factor = 0.2;

  akantu::Mesh mesh(spatial_dimension);

  akantu::MeshPartition * partition = NULL;
  if(prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("square_structured.msh", mesh);
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);

    // create the partition
    akantu::Vector<akantu::Int> part_tab(0,1);
    mesh.computeBoundingBox();
    akantu::Real rank_border = 0.5 * (mesh.getYMax() - mesh.getYMin());
    
    akantu::Mesh::type_iterator it  = mesh.firstType(spatial_dimension);
    akantu::Mesh::type_iterator end = mesh.lastType(spatial_dimension);
    for(; it != end; ++it) {
      akantu::ElementType c_type = *it;
      akantu::UInt nb_element = mesh.getNbElement(*it);
      
      for (akantu::UInt el=0; el<nb_element; ++el) {
	akantu::Real barycenter[spatial_dimension];
	mesh.getBarycenter(el,c_type,barycenter);
	if (barycenter[1] > rank_border)
	  part_tab.push_back(0);
	else
	  part_tab.push_back(1);
      }
    }

    partition->fillPartitionInformation(mesh,part_tab.storage());
  }

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);
  model->initParallel(partition);

  akantu::UInt nb_element = mesh.getNbElement(type);
  for(akantu::UInt i=0; i<nb_element; ++i) 
    proc_rank.push_back(prank);

  /// model initialization
  model->initVectors();
  
  /// set vectors to 0
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();

  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

  if(prank == 0)
    std::cout << model->getMaterial(0) << std::endl;

  model->setPBC(1,0,0);
  model->initPBC();
  model->assembleMassLumped();

  /// boundary conditions
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & coords = model->getFEM().getMesh().getNodes();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(std::abs(coords(i,1)-mesh.getYMax()) <= eps 
       || std::abs(coords(i,1)-mesh.getYMin()) <= eps) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  model->updateResidual();
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    if(s % 20 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
  }

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(iohelper::Dumper & dumper, const akantu::SolidMechanicsModel & model) {
  akantu::Real spatial_dimension = model.getSpatialDimension();
  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  akantu::UInt nb_element = model.getFEM().getMesh().getNbElement(type);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "pbc_parallel");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model.getMass().values,
			  spatial_dimension, "mass");
  dumper.AddNodeDataField(model.getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.AddElemDataField(proc_rank.values, 1, "rank");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(iohelper::Dumper & dumper) {
  dumper.Dump();
}
#endif
