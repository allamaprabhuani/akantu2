/**
 * @file   test_cohesive_extrinsic_IG_TG.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Thu April 18 13:31:00 2013
 *
 * @brief  Test for considering different cohesive properties for intergranular (IG) and
 * transgranular (TG) fractures in extrinsic cohesive elements
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
#if defined(AKANTU_USE_IOHELPER)
#  include "io_helper.hh"
#endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  const ElementType type = _triangle_6;

  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat", _explicit_dynamic, true, false);
  Real time_step = model.getStableTimeStep()*0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  //  std::cout << mesh << std::endl;

  const Mesh & mesh_facets = model.getMeshFacets();

  const ElementType type_facet = mesh.getFacetType(type);
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  Array<Real> & position = mesh.getNodes();

  Array<bool> & facet_check = model.getFacetsCheck();

  Real * bary_facet = new Real[spatial_dimension];
// first, the tag which shows grain ID should be read for each element

  // Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");
  // MeshUtils::buildAllFacets(mesh, mesh_facets);
  UInt nb_TG = 0;
  UInt nb_IG = 0;

  const Array< std::vector<Element> > & element_to_subelement = mesh_facets.getElementToSubelement(type_facet);
  Array<UInt> & facet_mat_by_type = model.getFacetMaterial(type_facet);
  ////////////////////////////////////////////

  for (UInt f = 0; f < nb_facet; ++f) {
    mesh_facets.getBarycenter(f, type_facet, bary_facet);
    if ((bary_facet[1] < 0.1 && bary_facet[1] > -0.1) ||(bary_facet[0] < 0.1 && bary_facet[0] > -0.1)) {
      facet_check(f) = true;
      const Element & el1 = element_to_subelement(f)[0];
      const Element & el2 = element_to_subelement(f)[1];
      UInt grain_id1 = mesh.getData<UInt>(el1.type, "tag_0")(el1.element);
      if(el2 != ElementNull) {
	UInt grain_id2 = mesh.getData<UInt>(el2.type, "tag_0")(el2.element);
	if (grain_id1 == grain_id2){
	  //transgranular = 0 indicator
	  facet_mat_by_type(f) = 1;
	  nb_TG++;
	} else  {
	  //intergranular = 1 indicator
	  facet_mat_by_type(f) = 2;
	  nb_IG++;
	}
      }
      std::cout << f << std::endl;
    }
    else
      facet_check(f) = false;
  }
  delete[] bary_facet;

  model.initFacetFilter();


   //  for ( UInt i = 0; i < nb_facet; ++i){

  // }
  //  std::cout << mesh << std::endl;

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99|| position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  model.updateResidual();

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();

  /// initial conditions
  Real loading_rate = 0.1;
  // bar_height  = 2
  Real VI = loading_rate * 2* 0.5;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
    velocity(n, 0) = loading_rate * position(n, 0);
  }

  // std::ofstream edis("edis.txt");
  // std::ofstream erev("erev.txt");

  //  Array<Real> & residual = model.getResidual();
  model.dump();
  //  const Array<Real> & stress = model.getMaterial(0).getStress(type);
  Real dispy = 0;
  // UInt nb_coh_elem = 0;

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    dispy += VI * time_step;
    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 1) > 0.99){
	displacement(n, 1) = dispy;
	velocity(n,1) = VI;}
      if (position(n, 1) < -0.99){
	displacement(n, 1) = -dispy;
	velocity(n,1) = -VI;}
      if (position(n, 0) > 0.99){
	displacement(n, 0) = dispy;
	velocity(n,0) = VI;}
      if (position(n, 0) < -0.99){
	displacement(n, 0) = -dispy;
	velocity(n,0) = -VI;}
    }

    model.checkCohesiveStress();

    model.explicitPred();
    // nb_coh_elem = mesh.getNbElement (FEM::getCohesiveElementType(type_facet));
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    model.dump();
    if(s % 10 == 0) {
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }

    // Real Ed = model.getEnergy("dissipated");

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;

  }

  // edis.close();
  // erev.close();

  //  mesh.write("mesh_final.msh");

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 40;

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  // for (UInt n = 0; n < position.getSize(); ++n) {
  //   for (UInt s = 0; s < spatial_dimension; ++s) {
  //     position(n, s) += displacement(n, s);
  //   }
  // }


  finalize();

  std::cout << "OK: test_cohesive_extrinsic_IG_TG was passed!" << std::endl;
  return EXIT_SUCCESS;
}
