/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
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
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

//#define CHECK_STRESS

int main(int argc, char *argv[])
{
#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_triangle_6;
  akantu::UInt paraview_type = TRIANGLE2;
#endif //AKANTU_USE_IOHELPER
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 10000;
  akantu::Real time_factor = 0.2;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar2.msh", mesh);

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


  model->readMaterials("material.dat");
  model->initMaterials();

  model->initModel();

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
    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] >= 9)
      model->getDisplacement().values[spatial_dimension*i] = (model->getFEM().getMesh().getNodes().values[spatial_dimension*i] - 9) / 100. ;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i] <= eps)
	model->getBoundary().values[spatial_dimension*i] = true;

    if(model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] <= eps ||
       model->getFEM().getMesh().getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  /// set the time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);


#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
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
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER


#ifdef CHECK_STRESS
  std::ofstream outfile;
  outfile.open("stress");
#endif // CHECK_STRESS

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef CHECK_STRESS
    /// search the position of the maximum of stress to determine the wave speed
    akantu::Real max_stress = std::numeric_limits<akantu::Real>::min();
    akantu::Real * stress = model->getMaterial(0).getStress(type).values;
    for (akantu::UInt i = 0; i < nb_element; ++i) {
      if(max_stress < stress[i*spatial_dimension*spatial_dimension]) {
	max_stress = stress[i*spatial_dimension*spatial_dimension];
      }
    }

    akantu::Real * coord    = model->getFEM().getMesh().getNodes().values;
    akantu::Real * disp_val = model->getDisplacement().values;
    akantu::UInt * conn     = model->getFEM().getMesh().getConnectivity(type).values;
    akantu::UInt nb_nodes_per_element = model->getFEM().getMesh().getNbNodesPerElement(type);
    akantu::Real * coords = new akantu::Real[spatial_dimension];
    akantu::Real min_x = std::numeric_limits<akantu::Real>::max();
    akantu::Real max_x = std::numeric_limits<akantu::Real>::min();

    akantu::Real stress_range = 5e7;
    for (akantu::UInt el = 0; el < nb_element; ++el) {
      if(stress[el*spatial_dimension*spatial_dimension] > max_stress - stress_range) {
	akantu::UInt el_offset  = el * nb_nodes_per_element;
	memset(coords, 0, spatial_dimension*sizeof(akantu::Real));
	for (akantu::UInt n = 0; n < nb_nodes_per_element; ++n) {
	  for (akantu::UInt i = 0; i < spatial_dimension; ++i) {
	    akantu::UInt node = conn[el_offset + n] * spatial_dimension;
	    coords[i] += (coord[node + i] + disp_val[node + i])
	      / ((akantu::Real) nb_nodes_per_element);
	  }
	}
	min_x = min_x < coords[0] ? min_x : coords[0];
	max_x = max_x > coords[0] ? max_x : coords[0];
      }
    }


    outfile << s << " " << .5 * (min_x + max_x) << " " << min_x << " " << max_x << " " << max_x - min_x << " " << max_stress << std::endl;

    delete [] coords;
#endif // CHECK_STRESS


#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

#ifdef CHECK_STRESS
  outfile.close();
#endif // CHECK_STRESS

  delete model;

  akantu::finalize();

  return EXIT_SUCCESS;
}
