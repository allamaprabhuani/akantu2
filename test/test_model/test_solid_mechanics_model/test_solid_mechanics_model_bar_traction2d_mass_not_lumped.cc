/**
 * @file   test_solid_mechanics_model.cc
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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"
using namespace iohelper;
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

//#define CHECK_STRESS
akantu::ElementType type = akantu::_triangle_3;
#ifdef AKANTU_USE_IOHELPER
  ElemType paraview_type = TRIANGLE1;
#endif //AKANTU_USE_IOHELPER

akantu::SolidMechanicsModel * model;
akantu::UInt spatial_dimension = 2;
akantu::UInt nb_nodes;
akantu::UInt nb_element;

akantu::Vector<akantu::Real> * lumped;

#ifdef AKANTU_USE_IOHELPER
static void paraviewInit(Dumper & dumper);
static void paraviewDump(Dumper & dumper);
#endif

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  //  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("bar1.msh", mesh);

  model = new akantu::SolidMechanicsModel(mesh);

  nb_nodes = model->getFEM().getMesh().getNbNodes();
  nb_element = model->getFEM().getMesh().getNbElement(type);

  lumped = new akantu::Vector<akantu::Real>(nb_nodes, spatial_dimension);

  /// model initialization
  model->initVectors();

  /// set vectors to 0
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();

  model->initExplicit();
  //model->initImplicit(true);
  model->initModel();
  model->readMaterials("material.dat");

  std::cout << model->getMaterial(0) << std::endl;

  model->initMaterials();

  model->initSolver();

  model->assembleMass();
  //  model->assembleStiffnessMatrix();

  model->getMassMatrix().lump(*lumped);

  /// boundary conditions
  akantu::Real eps = 1e-16;
  const akantu::Vector<akantu::Real> & pos = mesh.getNodes();
  akantu::Vector<akantu::Real> & disp = model->getDisplacement();
  akantu::Vector<bool> & boun = model->getBoundary();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  /// set the time step
  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  model->updateResidual();
  model->initialAcceleration();


#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  paraviewInit(dumper);
#endif //AKANTU_USE_IOHELPER

  std::ofstream energy;
  energy.open("energy_bar_2d_not_lumped.csv");
  energy << "id,rtime,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    //    model->implicitPred();
    // /// convergence loop
    // UInt count = 0;
    // Real error = 0.;
    // do {
    // 	std::cout << "passing step " << s << " " << s * time_step << "s - " << std::setw(4) << count << " : " << std::scientific << error << "\r" << std::flush;
    //   model->updateResidual();
    //   model->solveDynamic();
    //   model->implicitCorr();
    //   count++;
    // } while(!model->testConvergenceIncrement(1e-12, error) && (count < 1000));

    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    akantu::Real epot = model->getPotentialEnergy();
    akantu::Real ekin = model->getKineticEnergy();

    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 1 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER
    if(s % 100 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

  delete model;

  akantu::finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* Dumper vars                                                                */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER
void paraviewInit(Dumper & dumper) {
  dumper.SetMode(TEXT);
  dumper.SetPoints(model->getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, "bar2d_mass_not_lumped");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(type).values,
			 paraview_type, nb_element, C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values,
			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values,
			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model->getAcceleration().values,
			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model->getResidual().values,
			  spatial_dimension, "force");
  dumper.AddNodeDataField(model->getForce().values,
			  spatial_dimension, "applied_force");
  dumper.AddNodeDataField(lumped->values,
			  spatial_dimension, "mass");


  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
void paraviewDump(Dumper & dumper) {
  dumper.Dump();
}
#endif
