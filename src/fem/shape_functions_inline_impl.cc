/**
 * @file   shape_functions_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date   Fri Jul 15 19:41:58 2011
 *
 * @brief  ShapeFunctions inline implementation
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

__END_AKANTU__
#include "fem.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline UInt ShapeFunctions::getShapeSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_size = 0;
#define GET_SHAPE_SIZE(type)				\
  shape_size = ElementClass<type>::getShapeSize()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_SIZE);// ,

#undef GET_SHAPE_SIZE

  AKANTU_DEBUG_OUT();
  return shape_size;
}

/* -------------------------------------------------------------------------- */
inline UInt ShapeFunctions::getShapeDerivativesSize(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt shape_derivatives_size = 0;
#define GET_SHAPE_DERIVATIVES_SIZE(type)				\
  shape_derivatives_size = ElementClass<type>::getShapeDerivativesSize()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_SIZE);// ,

#undef GET_SHAPE_DERIVATIVES_SIZE

  AKANTU_DEBUG_OUT();
  return shape_derivatives_size;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::setControlPointsByType(const types::Matrix<Real> & points,
					    const GhostType & ghost_type) {
  control_points(type, ghost_type).shallowCopy(points);
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::interpolateElementalFieldOnControlPoints(const Vector<Real> &u_el,
							      Vector<Real> &uq,
							      GhostType ghost_type,
							      const Vector<Real> & shapes,
							      const Vector<UInt> * filter_elements) const {
  UInt nb_element;
  UInt nb_points = control_points(type, ghost_type).cols();
  UInt nb_nodes_per_element = ElementClass<type>::getShapeSize();
  UInt nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  Vector<Real>::const_iterator< types::Matrix<Real> > N_it;
  Vector<Real>::const_iterator< types::Matrix<Real> > u_it;
  Vector<Real>::iterator< types::Matrix<Real> > inter_u_it;

  Vector<Real> * filtered_N = NULL;
  if(filter_elements) {
    nb_element = filter_elements->getSize();
    filtered_N = new Vector<Real>(0, shapes.getNbComponent());
    FEM::filterQuadraturePointsData(*mesh, shapes, *filtered_N, type, ghost_type, filter_elements);
    N_it = filtered_N->begin_reinterpret(nb_nodes_per_element, nb_points, nb_element);
  } else {
    nb_element = mesh->getNbElement(type,ghost_type);
    N_it = shapes.begin_reinterpret(nb_nodes_per_element, nb_points, nb_element);
  }

  uq.resize(nb_element*nb_points);

  u_it       = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  inter_u_it = uq.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++N_it, ++u_it, ++inter_u_it) {
    const types::Matrix<Real> & u = *u_it;
    const types::Matrix<Real> & N = *N_it;
    types::Matrix<Real> & inter_u = *inter_u_it;

    inter_u.mul<false, false>(u, N);
  }

  if(filtered_N) delete filtered_N;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void ShapeFunctions::gradientElementalFieldOnControlPoints(const Vector<Real> &u_el,
							   Vector<Real> &out_nablauq,
							   GhostType ghost_type,
							   const Vector<Real> & shapes_derivatives,
							   const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element  = ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt nb_points             = control_points(type, ghost_type).cols();
  UInt element_dimension     = ElementClass<type>::getNaturalSpaceDimension();
  UInt nb_degree_of_freedom = u_el.getNbComponent() / nb_nodes_per_element;

  Vector<Real>::const_iterator< types::Matrix<Real> > B_it;
  Vector<Real>::const_iterator< types::Matrix<Real> > u_it;
  Vector<Real>::iterator< types::Matrix<Real> > nabla_u_it;

  UInt nb_element;
  Vector<Real> * filtered_B = NULL;
  if(filter_elements) {
    nb_element = filter_elements->getSize();
    filtered_B = new Vector<Real>(0, shapes_derivatives.getNbComponent());
    FEM::filterQuadraturePointsData(*mesh, shapes_derivatives, *filtered_B, type, ghost_type, filter_elements);
    B_it = filtered_B->begin(element_dimension, nb_nodes_per_element);
  } else {
    B_it = shapes_derivatives.begin(element_dimension, nb_nodes_per_element);
    nb_element = mesh->getNbElement(type, ghost_type);
  }

  out_nablauq.resize(nb_element*nb_points);

  u_it = u_el.begin(nb_degree_of_freedom, nb_nodes_per_element);
  nabla_u_it = out_nablauq.begin(nb_degree_of_freedom, element_dimension);

  for (UInt el = 0; el < nb_element; ++el, ++u_it) {
    const types::Matrix<Real> & u = *u_it;
    for (UInt q = 0; q < nb_points; ++q, ++B_it, ++nabla_u_it) {
      const types::Matrix<Real> & B = *B_it;
      types::Matrix<Real> & nabla_u = *nabla_u_it;

      nabla_u.mul<false, true>(u, B);
    }
  }

  if(filtered_B) delete filtered_B;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

