/**
 * @file   test_material_elastic_caughey.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed May  4 16:22:17 2011
 *
 * @brief  test of the material elastic caughey
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

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

int main(int argc, char *argv[])
{
#ifdef AKANTU_USE_IOHELPER
  akantu::ElementType type = akantu::_triangle_3;
#endif //AKANTU_USE_IOHELPER
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::initialize(&argc, &argv);

  akantu::debug::setDebugLevel(akantu::dblWarning);
  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);

  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  akantu::UInt nb_nodes = mesh.getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  akantu::UInt nb_element = mesh.getNbElement(type);
#endif //AKANTU_USE_IOHELPER
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


  model->initModel();
  model->readMaterials("material_elastic_caughey.dat");
  model->initMaterials();

  std::cout << model->getMaterial(0) << std::endl;

  model->assembleMassLumped();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    // model->getDisplacement().values[spatial_dimension*i]
    //   = mesh.getNodes().values[spatial_dimension*i] / 100.;
    if(mesh.getNodes().values[spatial_dimension*i] >= 9.9)
      model->getDisplacement().values[spatial_dimension*i]
     	= (mesh.getNodes().values[spatial_dimension*i] - 9) / 100. / 5.;

    if(mesh.getNodes().values[spatial_dimension*i] <= eps)
      model->getBoundary().values[spatial_dimension*i] = true;

    if(mesh.getNodes().values[spatial_dimension*i + 1] <= eps ||
       mesh.getNodes().values[spatial_dimension*i + 1] >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  model->synchronizeBoundaries();


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(mesh.getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d_mat_caughey");
  dumper.SetConnectivity((int *)mesh.getConnectivity(type).values,
			 TRIANGLE1, nb_element, C_MODE);
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
  dumper.AddElemDataField(dynamic_cast<akantu::MaterialElasticCaughey & >(model->getMaterial(0)).getStressViscosity(type).values,
			  spatial_dimension*spatial_dimension, "stress_visc");
  dumper.AddElemDataField(dynamic_cast<akantu::MaterialElasticCaughey & >(model->getMaterial(0)).getStressElastic(type).values,
			  spatial_dimension*spatial_dimension, "stress_el");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  /// Setting time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(akantu::UInt s = 1; s <= max_steps; ++s) {

    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    if(s % 10 == 0) {
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

  delete model;

  //  akantu::finalize();
  return EXIT_SUCCESS;
}
