/**
 * @file   test_cohesive_intrinsic_quadrangle.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Oct 03 10:20:53 2012
 *
 * @brief  Intrinsic cohesive elements' test for quadrangles
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

  const UInt spatial_dimension = 2;
  const UInt max_steps = 350;

  const ElementType type = _quadrangle_8;

  Mesh mesh(spatial_dimension);
  mesh.read("quadrangle.msh");


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");
  MeshUtils::buildAllFacets(mesh, mesh_facets);

  debug::setDebugLevel(dblDump);
  std::cout << mesh << std::endl;
  std::cout << mesh_facets << std::endl;
  debug::setDebugLevel(dblWarning);

  const ElementType type_facet = Mesh::getFacetType(type);

  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  //  const Vector<Real> & position = mesh.getNodes();
  //  Vector<Real> & displacement = model.getDisplacement();
  //  const Vector<UInt> & connectivity = mesh_facets.getConnectivity(type_facet);

  Vector<UInt> facet_insertion;
  Real * bary_facet = new Real[spatial_dimension];
  for (UInt f = 0; f < nb_facet; ++f) {
    mesh_facets.getBarycenter(f, type_facet, bary_facet);
    if (bary_facet[0] > -0.01 && bary_facet[0] < 0.01) facet_insertion.push_back(f);
  }
  delete[] bary_facet;

  MeshUtils::insertIntrinsicCohesiveElements(mesh,
					     mesh_facets,
					     type_facet,
					     facet_insertion);

  mesh.write("mesh_cohesive.msh");

  //  std::cout << mesh << std::endl;

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */


  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat");
  Real time_step = model.getStableTimeStep()*0.8;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


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

  model.setBaseName("intrinsic_quadrangle");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("force");
  model.dump();

  /// update displacement
  Vector<UInt> elements;
  Real * bary = new Real[spatial_dimension];
  for (UInt el = 0; el < nb_element; ++el) {
    mesh.getBarycenter(el, type, bary);
    if (bary[0] > 0.) elements.push_back(el);
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
      model.dump();
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

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 2;

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    return EXIT_FAILURE;
  }

  finalize();

  std::cout << "OK: test_cohesive_intrinsic_quadrangle was passed!" << std::endl;
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
  update.clear();

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
