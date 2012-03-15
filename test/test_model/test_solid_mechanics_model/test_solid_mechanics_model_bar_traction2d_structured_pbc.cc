/**
 * @file   test_solid_mechanics_model_bar_traction2d_structured_pbc.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Thu Dec 22 13:37:06 2011
 *
 * @brief  test of pbc class for SolidMechanicsModel
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

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  akantu::initialize(argc, argv);
  akantu::UInt spatial_dimension = 2;
  akantu::UInt max_steps = 1000;
  akantu::Real time_factor = 0.2;

  akantu::Real epot, ekin;

  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("bar_structured1.msh", mesh);

  akantu::SolidMechanicsModel * model = new akantu::SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  
  /// set vectors to 0
  akantu::UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getVelocity().values,     0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getAcceleration().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));
  memset(model->getDisplacement().values, 0,
	 spatial_dimension*nb_nodes*sizeof(akantu::Real));

  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

  std::cout << model->getMaterial(0) << std::endl;

  model->setPBC(1,0,0);
  model->initPBC();
  model->assembleMassLumped();

  /// boundary conditions
  mesh.computeBoundingBox();
  akantu::Real eps = 1e-16;
  akantu::Real signal_start = 0.6*mesh.getXMax();
  akantu::Real signal_end = 0.7*mesh.getXMax();
  akantu::Real delta_d = signal_end - signal_start;
  akantu::Real signal = 1.;
  const akantu::Vector<akantu::Real> & coords = model->getFEM().getMesh().getNodes();
  akantu::Vector<akantu::Real> & disp = model->getDisplacement();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(coords(i,0) >= signal_start && coords(i,0) <= signal_end) {
      if (coords(i,0) <= 0.5 * (signal_start + signal_end))
	disp(i,0) = (coords(i,0) - signal_start) * 2 * signal / delta_d;
      else
	disp(i,0) = (signal_end - coords(i,0)) * 2 * signal / delta_d;
    }
    
    if(coords(i,1) <= eps || coords(i,1) >= 1 - eps ) {
      model->getBoundary().values[spatial_dimension*i + 1] = true;
    }
  }

  akantu::Real time_step = model->getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model->setTimeStep(time_step);

  std::ofstream energy;
  energy.open("energy_2d_pbc.csv");
  energy << "id,epot,ekin,tot" << std::endl;

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  model->updateResidual();
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, *model);
#endif //AKANTU_USE_IOHELPER

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model->explicitPred();
    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

#ifdef AKANTU_USE_IOHELPER
    if(s % 20 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;
  }

  energy.close();

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
		   spatial_dimension, nb_nodes, "bar_traction_2d_structured_pbc");
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
