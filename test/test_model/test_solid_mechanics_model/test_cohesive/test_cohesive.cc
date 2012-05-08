/**
 * @file   test_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Feb 24 14:32:31 2012
 *
 * @brief  Test for cohesive elements
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
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 2500;

  const ElementType type = _triangle_3;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("mesh.msh", mesh);

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initExtrinsic("material.dat");
  Real time_step = model.getStableTimeStep()*0.8;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  const ElementType type_facet = mesh.getFacetElementType(type);

  const Mesh & mesh_facets = model.getMeshFacets();

  Vector<bool> & boundary = model.getBoundary();
  //  Vector<Real> & velocity     = model.getVelocity();
  Vector<Real> & displacement = model.getDisplacement();
  const Vector<Real> & position = mesh.getNodes();
  // Vector<Real> & acceleration = model.getAcceleration();
  // Vector<Real> & increment_acceleration = model.getIncrementAcceleration();
  // Vector<Real> & force        = model.getResidual();
  // const Vector<Real> & strain = model.getMaterial(0).getStrain(type);
  // const Vector<Real> & stress = model.getMaterial(0).getStress(type);

  UInt nb_nodes = model.getFEM().getMesh().getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 0) == -1 || position(n, 0) == 1) {
      boundary(n, 0) = boundary(n, 1) = true;
    }
  }


  iohelper::ElemType paraview_type = iohelper::TRIANGLE1;
  UInt nb_element = model.getFEM().getMesh().getNbElement(type);

  model.updateResidual();

  /// initialize the paraview output
  iohelper::DumperParaview dumper;
  //  dumper.SetMode(iohelper::TEXT);
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
  		   spatial_dimension, mesh.getNbNodes(), "explicit");
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
  			 paraview_type, nb_element, iohelper::C_MODE);
  dumper.AddNodeDataField(model.getDisplacement().values,
  			  spatial_dimension, "displacements");
  dumper.AddNodeDataField(model.getVelocity().values,
  			  spatial_dimension, "velocity");
  dumper.AddNodeDataField(model.getAcceleration().values,
  			  spatial_dimension, "acceleration");
  dumper.AddNodeDataField(model.getForce().values,
  			  spatial_dimension, "applied_force");
  dumper.AddNodeDataField(model.getResidual().values,
  			  spatial_dimension, "forces");
  dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
			  spatial_dimension*spatial_dimension, "stress");
  dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
  			  spatial_dimension*spatial_dimension, "strain");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("applied_force", 1);
  dumper.SetEmbeddedValue("forces", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();


  /// update displacement

  Real increment = 0.0001;

  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 0) == 1)
      displacement(n, 0) += increment;
  }


  std::ofstream edis("edis.txt");
  std::ofstream erev("erev.txt");

  /// main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.checkCohesiveStress();

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    if(s % 1 == 0) {

      dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		       spatial_dimension, mesh.getNbNodes(), "explicit");
      dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
			     paraview_type, nb_element, iohelper::C_MODE);
      dumper.AddNodeDataField(model.getDisplacement().values,
			      spatial_dimension, "displacements");
      dumper.AddNodeDataField(model.getVelocity().values,
			      spatial_dimension, "velocity");
      dumper.AddNodeDataField(model.getAcceleration().values,
			      spatial_dimension, "acceleration");
      dumper.AddNodeDataField(model.getForce().values,
			      spatial_dimension, "applied_force");
      dumper.AddNodeDataField(model.getResidual().values,
			      spatial_dimension, "forces");
      dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
      			      spatial_dimension*spatial_dimension, "strain");
      dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
       			      spatial_dimension*spatial_dimension, "stress");
      dumper.SetEmbeddedValue("displacements", 1);
      dumper.SetEmbeddedValue("applied_force", 1);
      dumper.SetEmbeddedValue("forces", 1);
      dumper.SetPrefix("paraview/");
      dumper.Init();
      dumper.Dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    /// update displacement
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (position(n, 0) == 1)
    	displacement(n, 0) += increment;
    }

    Real Ed = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getDissipatedEnergy();
    Real Er = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getReversibleEnergy();

    edis << s << " "
    	 << Ed << std::endl;

    erev << s << " "
    	 << Er << std::endl;

  }

  edis.close();
  erev.close();

  finalize();
  return EXIT_SUCCESS;
}
