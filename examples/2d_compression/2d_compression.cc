/**
 * @file   2d_compression.cc
 * @author Leo
 * @date   Thu Sep 30 17:05:00 2010
 *
 * @brief  2d dynamic explicit compression test with akantu
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

//using namespace akantu;

static void setInitialConditions(akantu::SolidMechanicsModel *);
static void setBoundaryConditions(akantu::SolidMechanicsModel *, akantu::Vector<akantu::UInt> *);
static void Apply_Displacement(akantu::SolidMechanicsModel *, akantu::Vector<akantu::UInt> *, akantu::Real);
static void reduceVelocities(akantu::SolidMechanicsModel *, akantu::Real);

#ifdef AKANTU_USE_IOHELPER
static void InitParaview(akantu::SolidMechanicsModel *);

DumperParaview dumper;
#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 30000;
  akantu::Real time_factor = 0.2;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("square_2d.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_elements = model->getFEM().getMesh().getNbElement(akantu::_triangle_3);

  /// model initialization
  model->initVectors();

  model->readMaterials("material.dat");
  model->initMaterials();

  model->initModel();
  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();


  /// Paraview Helper
  memset(model->getResidual().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getMaterial(0).getStrain(akantu::_triangle_3).values, 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(akantu::Real));
  memset(model->getMaterial(0).getStress(akantu::_triangle_3).values, 0,
	 spatial_dimension*spatial_dimension*nb_elements*sizeof(akantu::Real));


  /// boundary and initial conditions
  setInitialConditions(model);
  akantu::Vector<akantu::UInt> * face_node = new akantu::Vector<akantu::UInt>(0, 2, "face_node");
  setBoundaryConditions(model, face_node);

#ifdef AKANTU_USE_IOHELPER
  /// Paraview Helper
  InitParaview(model);
#endif //AKANTU_USE_IOHELPER

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    if(s<=100)
      Apply_Displacement(model, face_node, -0.01/100);

    model->explicitPred();

    model->updateResidual();

    model->updateAcceleration();
    model->explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    if(s % 200 == 0)
      dumper.Dump();
#endif //AKANTU_USE_IOHELPER

    if(s%100 == 0 && s>499)
      reduceVelocities(model, 0.95);

    if(s % 100 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  delete face_node;
  return EXIT_SUCCESS;
}



/**************************
 **   Static Functions   **
 **************************/

static void setInitialConditions(akantu::SolidMechanicsModel * model)
{
  akantu::UInt spatial_dimension = model->getSpatialDimension();
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  /// set vectors to 0
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
}


static void setBoundaryConditions(akantu::SolidMechanicsModel * model, akantu::Vector<akantu::UInt> * face_node)
{
  akantu::Real eps = 1e-16;
  akantu::Real * coord = model->getFEM().getMesh().getNodes().values;
  bool * id = model->getBoundary().values;
  akantu::Real y_min = std::numeric_limits<akantu::Real>::max();
  akantu::Real y_max = std::numeric_limits<akantu::Real>::min();
  akantu::UInt temp[2] = {0,0}, i;
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt spatial_dimension = model->getSpatialDimension();

  for(i = 0; i < nb_nodes; ++i) {
    y_min = (coord[spatial_dimension*i+1] < y_min) ? coord[spatial_dimension*i+1] : y_min;
    y_max = (coord[spatial_dimension*i+1] > y_max) ? coord[spatial_dimension*i+1] : y_max;
  }

  for (i = 0; i < nb_nodes; ++i) {
    if(coord[spatial_dimension*i+1]-eps <= y_min) {
      id[spatial_dimension*i+1] = true;
      temp[0] = 1; // face id
      temp[1] = i; // node id
      face_node->push_back(temp);
      continue;
    }

    if(coord[spatial_dimension*i+1]+eps >= y_max)
      id[spatial_dimension*i+1] = true;
  }

}

static void Apply_Displacement(akantu::SolidMechanicsModel * model, akantu::Vector<akantu::UInt> * face_node, akantu::Real delta)
{
  akantu::Real * disp_val = model->getDisplacement().values;
  const akantu::UInt nb_face_nodes = face_node->getSize();

  for(akantu::UInt i=0; i<nb_face_nodes; ++i)
    if(face_node->values[2*i] == 1)  { /// Node on top surface
      disp_val[2*(face_node->values[2*i+1])+1] += delta;
    }
}


// Artificial damping of velocities in order to reach a global static equilibrium
static void reduceVelocities(akantu::SolidMechanicsModel * model, akantu::Real ratio)
{
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::Real * velocities = model->getVelocity().values;

  if(ratio>1.) {
    fprintf(stderr,"**error** in Reduce_Velocities ratio bigger than 1!\n");
    exit(-1);
  }

  for(akantu::UInt i =0; i<nb_nodes; i++) {
    velocities[2*i] *= ratio;
    velocities[2*i+1] *= ratio;
  }
}

#ifdef AKANTU_USE_IOHELPER
static void InitParaview(akantu::SolidMechanicsModel * model)
{
  akantu::UInt spatial_dimension = model->getSpatialDimension();
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  akantu::UInt nb_elements = model->getFEM().getMesh().getNbElement(akantu::_triangle_3);

  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(akantu::_triangle_3).values,
			 TRIANGLE1, nb_elements, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(akantu::_triangle_3).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(akantu::_triangle_3).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}
#endif //AKANTU_USE_IOHELPER
