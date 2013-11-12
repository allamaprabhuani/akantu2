/**
 * @file   test_cohesive_extrinsic_tetrahedron.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Sep 12 11:50:14 2013
 *
 * @brief  Test for serial extrinsic cohesive elements for tetrahedron
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
#include "material_cohesive_linear.hh"
// #if defined(AKANTU_USE_IOHELPER)
// #  include "io_helper.hh"
// #endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real function(Real constant, Real x, Real y, Real z) {
  return constant + 2. * x + 3. * y + 4 * z;
}

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  // const UInt max_steps = 1000;
  // Real increment = 0.005;
  const UInt spatial_dimension = 3;

  ElementType type = _tetrahedron_10;
  ElementType type_facet = Mesh::getFacetType(type);
  ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);

  Mesh mesh(spatial_dimension);
  mesh.read("tetrahedron.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat", _explicit_lumped_mass, true);

  const MaterialCohesiveLinear<3> & mat_cohesive
    = dynamic_cast < const MaterialCohesiveLinear<3> & > (model.getMaterial(1));

  const Real sigma_c = mat_cohesive.getParam<Real>("sigma_c");
  const Real beta = mat_cohesive.getParam<Real>("beta");
  //  const Real G_cI = mat_cohesive.getParam<Real>("G_cI");

  Array<Real> & position = mesh.getNodes();

  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  /// compute quadrature points positions on facets
  const Mesh & mesh_facets = model.getMeshFacets();
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  UInt nb_quad_per_facet = model.getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

  Array<Real> quad_facets(nb_tot_quad, spatial_dimension);

  model.getFEM("FacetsFEM").interpolateOnQuadraturePoints(position,
							  quad_facets,
							  spatial_dimension,
							  type_facet);

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */


  /// compute quadrature points position of the elements
  UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_tot_quad_el = nb_quad_per_element * nb_element;

  Array<Real> quad_elements(nb_tot_quad_el, spatial_dimension);


  model.getFEM().interpolateOnQuadraturePoints(position,
					       quad_elements,
					       spatial_dimension,
					       type);

  /// assign some values to stresses
  Array<Real> & stress
    = const_cast<Array<Real>&>(model.getMaterial(0).getStress(type));

  Array<Real>::iterator<Matrix<Real> > stress_it
    = stress.begin(spatial_dimension, spatial_dimension);

  for (UInt q = 0; q < nb_tot_quad_el; ++q, ++stress_it) {

    /// compute values
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = i; j < spatial_dimension; ++j) {
    	UInt index = i * spatial_dimension + j;
    	(*stress_it)(i, j) = index * function(sigma_c * 5,
    					      quad_elements(q, 0),
    					      quad_elements(q, 1),
    					      quad_elements(q, 2));
      }
    }

    /// fill symmetrical part
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < i; ++j) {
    	(*stress_it)(i, j) = (*stress_it)(j, i);
      }
    }

    // stress_it->clear();
    // for (UInt i = 0; i < spatial_dimension; ++i)
    //   (*stress_it)(i, i) = sigma_c * 5;
  }


  /// compute stress on facet quads
  Array<Real> stress_facets(nb_tot_quad, spatial_dimension * spatial_dimension);

  Array<Real>::iterator<Matrix<Real> > stress_facets_it
    = stress_facets.begin(spatial_dimension, spatial_dimension);

  for (UInt q = 0; q < nb_tot_quad; ++q, ++stress_facets_it) {
    /// compute values
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = i; j < spatial_dimension; ++j) {
    	UInt index = i * spatial_dimension + j;
    	(*stress_facets_it)(i, j) = index * function(sigma_c * 5,
    						     quad_facets(q, 0),
    						     quad_facets(q, 1),
    						     quad_facets(q, 2));
      }
    }

    /// fill symmetrical part
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < i; ++j) {
    	(*stress_facets_it)(i, j) = (*stress_facets_it)(j, i);
      }
    }


    // stress_facets_it->clear();
    // for (UInt i = 0; i < spatial_dimension; ++i)
    //   (*stress_facets_it)(i, i) = sigma_c * 5;
  }

  // /// compute facet area
  // Array<Real> integration_constant(nb_tot_quad);
  // integration_constant.set(1.);

  // Real facets_area = model.getFEM("FacetsFEM").integrate(integration_constant,
  // 							 type_facet);

  // facets_area -= 2 * 2 * 6;


  /// insert cohesive elements
  model.checkCohesiveStress();


  /// check insertion stress
  const Array<Real> & normals =
    model.getFEM("FacetsFEM").getNormalsOnQuadPoints(type_facet);
  const Array<Real> & tangents = model.getTangents(type_facet);
  const Array<Real> & sigma_c_eff = mat_cohesive.getInsertionTraction(type_cohesive);

  Vector<Real> normal_stress(spatial_dimension);

  const Array<std::vector<Element> > & coh_element_to_facet
    = mesh_facets.getElementToSubelement(type_facet);

  Array<Real>::iterator<Matrix<Real> > quad_facet_stress
    = stress_facets.begin(spatial_dimension, spatial_dimension);

  Array<Real>::const_iterator<Vector<Real> > quad_normal
    = normals.begin(spatial_dimension);

  Array<Real>::const_iterator<Vector<Real> > quad_tangents
    = tangents.begin(tangents.getNbComponent());

  for (UInt f = 0; f < nb_facet; ++f) {
    const Element & cohesive_element = coh_element_to_facet(f)[1];

    for (UInt q = 0; q < nb_quad_per_facet; ++q, ++quad_facet_stress,
	   ++quad_normal, ++quad_tangents) {
      if (cohesive_element == ElementNull) continue;

      normal_stress.mul<false>(*quad_facet_stress, *quad_normal);

      Real normal_contrib = normal_stress.dot(*quad_normal);

      Real first_tangent_contrib = 0;

      for (UInt dim = 0; dim < spatial_dimension; ++dim)
	first_tangent_contrib += normal_stress(dim) * (*quad_tangents)(dim);

      Real second_tangent_contrib = 0;

      for (UInt dim = 0; dim < spatial_dimension; ++dim)
	second_tangent_contrib
	  += normal_stress(dim) * (*quad_tangents)(dim + spatial_dimension);

      Real tangent_contrib = std::sqrt(first_tangent_contrib * first_tangent_contrib +
				       second_tangent_contrib * second_tangent_contrib);

      normal_contrib = std::max(0., normal_contrib);

      Real effective_norm = std::sqrt(normal_contrib * normal_contrib
				      + tangent_contrib * tangent_contrib / beta / beta);

      if (effective_norm < sigma_c) continue;

      if (!Math::are_float_equal(effective_norm,
				 sigma_c_eff(cohesive_element.element
					     * nb_quad_per_facet + q))) {
	std::cout << "Insertion tractions do not match" << std::endl;
	finalize();
	return EXIT_FAILURE;
      }
    }
  }

  // /* ------------------------------------------------------------------------ */
  // /* Check dissipated energy                                                  */
  // /* ------------------------------------------------------------------------ */


  // model.setBaseName("extrinsic_tetrahedron");
  // model.addDumpFieldVector("displacement");
  // model.addDumpField("stress");
  // model.dump();

  // DumperParaview dumper("cohesive_elements_tetrahedron_fragmentation");
  // dumper.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_cohesive);
  // DumperIOHelper::Field * cohesive_displacement =
  //   new DumperIOHelper::NodalField<Real>(model.getDisplacement());
  // cohesive_displacement->setPadding(3);
  // dumper.registerField("displacement", cohesive_displacement);
  // // dumper.registerField("damage", new DumperIOHelper::
  // // 		       HomogenizedField<Real,
  // // 					DumperIOHelper::InternalMaterialField>(model,
  // // 									       "damage",
  // // 									       spatial_dimension,
  // // 									       _not_ghost,
  // // 									       _ek_cohesive));

  // dumper.dump();

  // /// update displacement
  // UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  // Real * bary = new Real[spatial_dimension];

  // const Array<UInt> & connectivity = mesh.getConnectivity(type);
  // Array<Real> & displacement = model.getDisplacement();

  // for (UInt s = 1; s <= max_steps; ++s) {

  //   model.updateResidual();

  //   for (UInt el = 0; el < nb_element; ++el) {
  //     mesh.getBarycenter(el, type, bary);
  //     for (UInt n = 0; n < nb_nodes_per_element; ++n) {
  // 	UInt node = connectivity(el, n);

  // 	for (UInt dim = 0; dim < spatial_dimension; ++dim) {
  // 	  displacement(node, dim) += increment * bary[dim];
  // 	}
  //     }
  //   }

  //   if (s % 100 == 0) {
  //     model.dump();
  //     dumper.dump();
  //   }
  // }

  // delete[] bary;


  // Real theoretical_Ed = facets_area * G_cI;
  // Real Ed = model.getEnergy("dissipated");

  // std::cout << Ed << " " << theoretical_Ed << std::endl;

  // Array<Real> & velocity = model.getVelocity();
  // Array<bool> & boundary = model.getBoundary();
  // Array<Real> & displacement = model.getDisplacement();
  // //  const Array<Real> & residual = model.getResidual();

  // UInt nb_nodes = mesh.getNbNodes();

  // /// boundary conditions
  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
  //     boundary(n, 0) = true;
  // }

  // model.updateResidual();

  // model.setBaseName("extrinsic_tetrahedron");
  // model.addDumpFieldVector("displacement");
  // model.addDumpField("velocity"    );
  // model.addDumpField("acceleration");
  // model.addDumpField("residual"    );
  // model.addDumpField("stress");
  // model.addDumpField("strain");
  // model.dump();

  // /// initial conditions
  // Real loading_rate = 0.5;
  // Real disp_update = loading_rate * time_step;
  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   velocity(n, 0) = loading_rate * position(n, 0);
  // }

  // //  const Array<Real> & stress = model.getMaterial(0).getStress(type);

  // /// Main loop
  // for (UInt s = 1; s <= max_steps; ++s) {

  //   /// update displacement on extreme nodes
  //   for (UInt n = 0; n < nb_nodes; ++n) {
  //     if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
  // 	displacement(n, 0) += disp_update * position(n, 0);
  //   }

  //   model.checkCohesiveStress();

  //   model.explicitPred();
  //   model.updateResidual();
  //   model.updateAcceleration();
  //   model.explicitCorr();

  //   if(s % 10 == 0) {
  //     model.dump();

  //     std::cout << "passing step " << s << "/" << max_steps << std::endl;
  //   }

  // }

  // //  mesh.write("mesh_final.msh");

  // Real Ed = model.getEnergy("dissipated");
  // Real Edt = 400;

  // std::cout << Ed << " " << Edt << std::endl;

  // // if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
  // //   std::cout << "The dissipated energy is incorrect" << std::endl;
  // //   return EXIT_FAILURE;
  // // }


  finalize();

  std::cout << "OK: test_cohesive_extrinsic was passed!" << std::endl;
  return EXIT_SUCCESS;
}
