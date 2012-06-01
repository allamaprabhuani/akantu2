/**
 * @file   test_material_viscoelastic.cc
 * @author Vlad Yastrebov <vladislav.yastrebov@epfl.ch> 
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Feb 9 2012 
 *
 * @brief  test of the material viscoelastic 
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
#include <sstream>
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
#  include "io_helper_tools.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::debug::setDebugLevel(akantu::dblWarning);

  const ElementType element_type = TYPE;
#ifdef AKANTU_USE_IOHELPER
  iohelper::ElemType ioh_type = getIOHelperType(element_type);
#endif //AKANTU_USE_IOHELPER
  const UInt dim = ElementClass<TYPE>::getSpatialDimension();

  /// load mesh
  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  std::stringstream meshname_sstr;
  meshname_sstr << "test_material_viscoelastic_" << element_type << ".msh";
  mesh_io.read(meshname_sstr.str().c_str(), mesh);

  UInt max_steps = 5000;
  Real time_factor = 0.8;
  UInt nb_nodes = mesh.getNbNodes();

  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initVectors();
  model.getForce().clear();
  model.getVelocity().clear();
  model.getAcceleration().clear();
  model.getDisplacement().clear();
  model.updateResidual();

  model.initExplicit();
  model.initModel();
  model.readMaterials("material_viscoelastic.dat");
  model.initMaterials();

  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  Real eps = 1e-16;
  for (UInt i = 0; i < nb_nodes; ++i) {
    if(mesh.getNodes().values[dim*i] >= 9)
      model.getDisplacement().values[dim*i+1]
     	= (mesh.getNodes().values[dim*i] - 9) / 100.;

    if(mesh.getNodes().values[dim*i] <= eps)
      model.getBoundary().values[dim*i] = true;

/*
    if(mesh.getNodes().values[dim*i + 1] <= eps ||
       mesh.getNodes().values[dim*i + 1] >= 1 - eps ) {
      model.getBoundary().values[dim*i + 1] = true;
    }
*/
  }

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  // paraviewInit(dumper, model, element_type, "test_viscoelastic2");
  dumper.SetMode(iohelper::BASE64);
  Material & mat = model.getMaterial(0);

  std::stringstream sstr; sstr << "test_material_viscoelastic_" << TYPE;

  dumper.SetPoints(model.getFEM().getMesh().getNodes().values, dim, nb_nodes, sstr.str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 ioh_type, model.getFEM().getMesh().getNbElement(TYPE), iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model.getForce().values, 2, "force");
  dumper.AddNodeDataField(model.getMass().values, 1, "Mass");
  dumper.AddNodeDataField(model.getResidual().values, 2, "residual");
  dumper.AddElemDataField(mat.getStrain(TYPE).values, 4, "strain");
  dumper.AddElemDataField(mat.getStress(TYPE).values, 4, "stress");

//  Real * dam = mat.getHistoryIntegral(TYPE).values;
//  dumper.AddElemDataField(dam, 4, "damage");

  try{
    const Vector<Real> & history = mat.getVector("history_integral", TYPE);
    dumper.AddElemDataField(history.storage(), 4, "history");
  } catch (...) {};

  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetEmbeddedValue("residual", 1);
  dumper.SetEmbeddedValue("velocity", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

#endif //AKANTU_USE_IOHELPER

  std::stringstream filename_sstr;
  filename_sstr << "test_material_viscoelastic_" << element_type << ".out";
  std::ofstream energy;
  energy.open(filename_sstr.str().c_str());
  energy << "id epot ekin tot" << std::endl;

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    Real epot = model.getPotentialEnergy();
    Real ekin = model.getKineticEnergy();

    if(s % 10 == 0) {
       std::cerr << "passing step " << s << "/" << max_steps << std::endl;
       energy << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;
    }

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER
  }

  finalize();
  return EXIT_SUCCESS;
}
