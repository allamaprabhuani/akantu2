/**
 * @file   colone_weight.cc
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @date   Wed Mar 23 15:42:36 2011
 *
 * @brief
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
#include "aka_vector.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[]) {

  // chose if you use hexahedron elements
  bool use_hexa = false;

  std::stringstream mesh_file;
  std::stringstream output;
  std::stringstream energy_file;
  akantu::ElementType type;
#ifdef AKANTU_USE_IOHELPER
   akantu::UInt paraview_type;
#endif //AKANTU_USE_IOHELPER
   UInt vel_damping_interval;

  if (use_hexa) {
    type = akantu::_hexahedron_8;
    mesh_file << "colone_hexa.msh";
    output << "paraview/test_weight_hexa";
    energy_file << "energy_hexa.csv";
#ifdef AKANTU_USE_IOHELPER
    paraview_type = HEX1;
#endif //AKANTU_USE_IOHELPER
    vel_damping_interval =4;
  }
  else {
    type = akantu::_tetrahedron_4;
    mesh_file << "colone_tetra.msh";
    output << "paraview/test_weight_tetra";
    energy_file << "energy_tetra.csv";
#ifdef AKANTU_USE_IOHELPER
    paraview_type = TETRA1;
#endif //AKANTU_USE_IOHELPER
    vel_damping_interval = 8;
  }

  akantu::UInt spatial_dimension = 3;

  akantu::UInt max_steps = 2000;
  akantu::Real time_factor = 0.8;

  akantu::initialize(&argc, &argv);

  //  akantu::Real epot, ekin;
  akantu::Mesh mesh(spatial_dimension);

  akantu::MeshIOMSH mesh_io;
  mesh_io.read(mesh_file.str().c_str(), mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

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
  model->readMaterials("material_colone.dat");
  model->initMaterials();

  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();


#ifdef AKANTU_USE_IOHELPER
  /// set to 0 only for the first paraview dump
  //  memset(model->getResidual().values, 0,
  //	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  // memset(model->getMaterial(0).getStrain(type).values, 0,
  // 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
  // memset(model->getMaterial(0).getStress(type).values, 0,
  // 	 spatial_dimension*spatial_dimension*nb_element*sizeof(akantu::Real));
#endif //AKANTU_USE_IOHELPER


  /// boundary conditions
  const akantu::Vector<Real> & position = model->getFEM().getMesh().getNodes();
  akantu::Vector<bool> & boundary = model->getBoundary();
  akantu::Vector<Real> & force = model->getForce();
  const akantu::Vector<Real> & mass = model->getMass();

  akantu::Real z_min = position(0, 2);
  for (unsigned int i = 0; i < nb_nodes; ++i) {
    if(position(i, 2) < z_min)
      z_min = position(i, 2);
  }

  akantu::Real eps = 1e-13;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(fabs(position(i, 2) - z_min) <= eps)
      boundary(i,2) = true;
    else
      force(i,2) = -mass(i,0) * 9.81;
  }

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);
  //  model->setTimeStep(3.54379e-07);

  model->updateResidual();
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
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(type).values,
   			  spatial_dimension*spatial_dimension, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(type).values,
   			  spatial_dimension*spatial_dimension, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix(output.str().c_str());
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  double total_time = 0.;
  akantu::Vector<Real> & velocity = model->getVelocity();

  std::ofstream energy;
  energy.open(energy_file.str().c_str());
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if (s % vel_damping_interval == 0) {
      for (akantu::UInt i = 0; i < nb_nodes; ++i) {
	velocity(i, 0) *= 0.9;
	velocity(i, 1) *= 0.9;
	velocity(i, 2) *= 0.9;
      }
    }

#ifdef AKANTU_USE_IOHELPER
    if(s % 1 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
    if(s % 10 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  akantu::finalize();

  return EXIT_SUCCESS;
}
