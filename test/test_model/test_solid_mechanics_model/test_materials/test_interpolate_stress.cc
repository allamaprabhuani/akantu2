/**
 * @file   test_interpolate_stress.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Wed Jun  6 14:25:54 2012
 *
 * @brief  Test for the stress interpolation function
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
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

static Real plane(Real, Real);

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const ElementType type = _triangle_6;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);
  const ElementType type_facet = mesh.getFacetElementType(type);;

  MeshUtils::buildAllFacets(mesh, mesh);

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull("../material.dat");

  Vector<Real> & position = mesh.getNodes();
  UInt nb_facet = mesh.getNbElement(type_facet);
  UInt nb_element = mesh.getNbElement(type);

  /// compute quadrature points positions on facets
  typedef FEMTemplate<IntegratorGauss,ShapeLagrange> MyFEMType;

  model.registerFEMObject<MyFEMType>("FacetsFEM", mesh, spatial_dimension-1);
  model.getFEM("FacetsFEM").initShapeFunctions();

  UInt nb_quad_per_facet = model.getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

  Vector<Real> quad_facets(nb_tot_quad, spatial_dimension);

  model.getFEM("FacetsFEM").interpolateOnQuadraturePoints(position,
							  quad_facets,
							  spatial_dimension,
							  type_facet);

  Vector<Element> & facet_to_element = mesh.getSubelementToElement(type);
  UInt nb_facet_per_elem = facet_to_element.getNbComponent();

  Vector<Real> el_q_facet(nb_element * nb_facet_per_elem * nb_quad_per_facet,
			  spatial_dimension);

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt f = 0; f < nb_facet_per_elem; ++f) {
      UInt global_facet = facet_to_element(el, f).element;

      for (UInt q = 0; q < nb_quad_per_facet; ++q) {
	for (UInt s = 0; s < spatial_dimension; ++s) {
	  el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
		     + f * nb_quad_per_facet + q, s)
	    = quad_facets(global_facet * nb_quad_per_facet + q, s);
	}
      }
    }
  }


  /// compute quadrature points position of the elements
  UInt nb_quad_per_element = model.getFEM().getNbQuadraturePoints(type);
  UInt nb_tot_quad_el = nb_quad_per_element * nb_element;

  Vector<Real> quad_elements(nb_tot_quad_el, spatial_dimension);


  model.getFEM().interpolateOnQuadraturePoints(position,
					       quad_elements,
					       spatial_dimension,
					       type);

  /// assign some values to stresses
  Vector<Real> & stress
    = const_cast<Vector<Real>&>(model.getMaterial(0).getStress(type));

  for (UInt q = 0; q < nb_tot_quad_el; ++q) {
    for (UInt s = 0; s < spatial_dimension * spatial_dimension; ++s) {
      stress(q, s) = s * plane(quad_elements(q, 0), quad_elements(q, 1));
    }
  }

  /// interpolate stresses on facets' quadrature points
  model.getMaterial(0).initInterpolateElementalField();

  Vector<Real> interpolated_stress(nb_element * nb_facet_per_elem * nb_quad_per_facet,
				   stress.getNbComponent());
  model.getMaterial(0).interpolateStress(type,
					 el_q_facet,
					 interpolated_stress);

  /// check results
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt f = 0; f < nb_facet_per_elem; ++f) {

      for (UInt q = 0; q < nb_quad_per_facet; ++q) {
	for (UInt s = 0; s < spatial_dimension * spatial_dimension; ++s) {

	  Real x = el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
			      + f * nb_quad_per_facet + q, 0);
	  Real y = el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet
			      + f * nb_quad_per_facet + q, 1);

	  Real theoretical = s * plane(x, y);

	  Real numerical = interpolated_stress(el * nb_facet_per_elem * nb_quad_per_facet
					       + f * nb_quad_per_facet + q, s);

	  Real tolerance = 1000 * (theoretical + 1) * std::numeric_limits<Real>::epsilon();

	  if (std::abs(theoretical - numerical) > tolerance) {
	    std::cout << "Theoretical and numerical values aren't coincident!" << std::endl;
	    return EXIT_FAILURE;
	  }

	}
      }
    }
  }

  std::cout << "OK: Stress interpolation test passed." << std::endl;
  return EXIT_SUCCESS;
}

Real plane(Real x, Real y) {
  return 1. + 2. * x + 3. * y;
}
