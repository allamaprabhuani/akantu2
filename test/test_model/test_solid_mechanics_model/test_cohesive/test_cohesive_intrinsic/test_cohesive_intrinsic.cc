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
#include <iostream>

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

static void updateDisplacement(SolidMechanicsModelCohesive &,
			       Vector<UInt> &,
			       ElementType,
			       Real);

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  //  debug::setDebugLevel(dblDump);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 350;

  const ElementType type = _triangle_6;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("mesh.msh", mesh);

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat");
  Real time_step = model.getStableTimeStep()*0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  //  std::cout << mesh << std::endl;

  const Mesh & mesh_facets = model.getMeshFacets();

  const ElementType type_facet = mesh.getFacetElementType(type);
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  //  const Vector<Real> & position = mesh.getNodes();
  //  Vector<Real> & displacement = model.getDisplacement();
  //  const Vector<UInt> & connectivity = mesh_facets.getConnectivity(type_facet);

  Vector<UInt> facet_insertion;
  Real * bary_facet = new Real[spatial_dimension];
  for (UInt f = 0; f < nb_facet; ++f) {
    mesh_facets.getBarycenter(f, type_facet, bary_facet);
    if (bary_facet[0] == -0.25) facet_insertion.push_back(f);
  }
  delete[] bary_facet;

  model.insertCohesiveElements(facet_insertion);

  //  mesh_io.write("mesh_cohesive.msh", mesh);

  //  std::cout << mesh << std::endl;

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  Vector<bool> & boundary = model.getBoundary();
  //  const Vector<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_element = mesh.getNbElement(type);

  /// boundary conditions
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    for (UInt n = 0; n < nb_nodes; ++n) {
      boundary(n, dim) = true;
    }
  }

  model.updateResidual();

  // iohelper::ElemType paraview_type = iohelper::TRIANGLE2;

  // /// initialize the paraview output
  // iohelper::DumperParaview dumper;
  // dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
  // 		   spatial_dimension, nb_nodes, "explicit");
  // dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(type).values,
  // 			 paraview_type, nb_element, iohelper::C_MODE);
  // dumper.AddNodeDataField(model.getDisplacement().values,
  // 			  spatial_dimension, "displacements");
  // dumper.AddNodeDataField(model.getVelocity().values,
  // 			  spatial_dimension, "velocity");
  // dumper.AddNodeDataField(model.getAcceleration().values,
  // 			  spatial_dimension, "acceleration");
  // dumper.AddNodeDataField(model.getForce().values,
  // 			  spatial_dimension, "applied_force");
  // dumper.AddNodeDataField(model.getResidual().values,
  // 			  spatial_dimension, "forces");
  // dumper.AddElemDataField(model.getMaterial(0).getStrain(type).values,
  // 			  spatial_dimension*spatial_dimension, "strain");
  // dumper.AddElemDataField(model.getMaterial(0).getStress(type).values,
  // 			  spatial_dimension*spatial_dimension, "stress");
  // dumper.SetEmbeddedValue("displacements", 1);
  // dumper.SetEmbeddedValue("applied_force", 1);
  // dumper.SetEmbeddedValue("forces", 1);
  // dumper.SetPrefix("paraview/");
  // dumper.Init();
  // dumper.Dump();


  /// update displacement
  Vector<UInt> elements;
  Real * bary = new Real[spatial_dimension];
  for (UInt el = 0; el < nb_element; ++el) {
    mesh.getBarycenter(el, type, bary);
    if (bary[0] > -0.25) elements.push_back(el);
  }
  delete[] bary;

  Real increment = 0.01;

  updateDisplacement(model, elements, type, increment);

  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   if (position(n, 1) + displacement(n, 1) > 0) {
  //     if (position(n, 0) == 0) {
  // 	displacement(n, 1) -= 0.25;
  //     }
  //     if (position(n, 0) == 1) {
  // 	displacement(n, 1) += 0.25;
  //     }
  //   }
  // }


  // std::ofstream edis("edis.txt");
  // std::ofstream erev("erev.txt");

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    updateDisplacement(model, elements, type, increment);

    if(s % 1 == 0) {
      //      dumper.Dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // // update displacement
    // for (UInt n = 0; n < nb_nodes; ++n) {
    //   if (position(n, 1) + displacement(n, 1) > 0) {
    // 	displacement(n, 0) -= 0.01;
    //   }
    // }

    //    Real Ed = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getDissipatedEnergy();
    //    Real Er = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;

  }

  // edis.close();
  // erev.close();

  Real Ed = model.getDissipatedEnergy();

  Real Edt = 2 * sqrt(2);

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  std::cout << "OK: test_cohesive_intrinsic was passed!" << std::endl;
  return EXIT_SUCCESS;
}


static void updateDisplacement(SolidMechanicsModelCohesive & model,
			       Vector<UInt> & elements,
			       ElementType type,
			       Real increment) {

  Mesh & mesh = model.getFEM().getMesh();
  UInt nb_element = elements.getSize();
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  const Vector<UInt> & connectivity = mesh.getConnectivity(type);
  Vector<Real> & displacement = model.getDisplacement();
  Vector<bool> update(nb_nodes);

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = connectivity(elements(el), n);
      if (!update(node)) {
	displacement(node, 0) += increment;
	//	displacement(node, 1) += increment;
	update(node) = true;
      }
    }
  }
}
