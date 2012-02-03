/**
 * @file   shape_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Thu Feb  2 15:55:15 2012
 *
 * @brief  
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

#include "aka_common.hh"
#include "cohesive_element.hh"
#include "shape_cohesive.hh"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

template <class ShapeFunction>
template <ElementType type>
void ShapeCohesive<ShapeFunction>::precomputeShapesOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const ElementType sub_type = CohesiveElementSubElementType<type>::type;

  Real * natural_coords = this->control_points(sub_type, ghost_type).storage();
  UInt nb_points        = this->control_points(sub_type, ghost_type).getSize();

  UInt size_of_shapes    = CohesiveElement<type>::getShapeSize();

  UInt nb_element = this->mesh->getConnectivity(type, ghost_type).getSize();;

  Vector<Real> & shapes_tmp = sub_type_shape_function.shapes.alloc(nb_element*nb_points,
								   size_of_shapes,
								   type,
								   ghost_type);

  Real * shapes_val    = shapes_tmp.storage();

  for (UInt elem = 0; elem < nb_element; ++elem) {
    ElementClass<sub_type>::computeShapes(natural_coords,
				      nb_points,shapes_val);

    shapes_val += size_of_shapes*nb_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
template <ElementType type>
void ShapeCohesive<ShapeFunction>::precomputeShapeDerivativesOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const ElementType sub_type = CohesiveElementSubElementType<type>::type;

  Real * coord = this->mesh->getNodes().storage();
  UInt spatial_dimension = this->mesh->getSpatialDimension();

  UInt nb_nodes_per_element = CohesiveElement<type>::getNbNodesPerElement(type);
  UInt nb_nodes_per_sub_element = ElementClass<sub_type>::getNbNodesPerElement();
  UInt size_of_shapesd      = CohesiveElement<type>::getShapeDerivativesSize();
  UInt nb_points            = this->control_points(type, ghost_type).getSize();
  Real * natural_coords     = this->control_points(type, ghost_type).storage();


  UInt * elem_val = this->mesh->getConnectivity(type, ghost_type).storage();;
  UInt nb_element = this->mesh->getConnectivity(type, ghost_type).getSize();

  Vector<Real> & shapes_derivatives_tmp =
    sub_type_shape_function.shapes_derivatives.alloc(nb_element*nb_points,
						     size_of_shapesd,
						     type,
						     ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.storage();
  Real local_coord[spatial_dimension * nb_nodes_per_element];

  CohesiveReduceFunctionMean reduce_function;

  for (UInt elem = 0; elem < nb_element; ++elem) {

    UInt el_offset = elem * nb_nodes_per_element;

    for (UInt n = 0; n < nb_nodes_per_sub_element; ++n) {
      for (UInt d = 0; d < spatial_dimension; ++d) {
	Real u_plus  = coord[elem_val[el_offset + n] * spatial_dimension];
	Real u_minus = coord[elem_val[el_offset + n + nb_nodes_per_sub_element] * spatial_dimension];
	local_coord[n*spatial_dimension + d] = reduce_function(u_plus, u_minus);
      }
    }


    sub_type_shape_function.
      computeShapeDerivativesOnCPointsByElement<sub_type>(spatial_dimension,
							  local_coord,
							  nb_nodes_per_element,
							  natural_coords,
							  nb_points,
							  shapesd_val);

    shapesd_val += size_of_shapesd*nb_points;
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
template <ElementType type, class ReduceFunction>
void ShapeCohesive<ShapeFunction>::interpolateOnControlPoints(const Vector<Real> &in_u,
					       Vector<Real> &out_uq,
					       UInt nb_degree_of_freedom,
					       GhostType ghost_type,
					       const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  const ElementType sub_type = CohesiveElementSubElementType<type>::type;

  AKANTU_DEBUG_ASSERT(sub_type_shape_function.shapes.exists(sub_type, ghost_type),
		      "No shapes for the type "
		      << sub_type_shape_function.shapes.printType(sub_type, ghost_type));

  UInt nb_element = this->mesh->getNbElement(type, ghost_type);
  UInt * conn_val = this->mesh->getConnectivity(type, ghost_type).storage();

  const Vector<Real> & shapes = sub_type_shape_function.shapes(sub_type, ghost_type);

  UInt nb_nodes_per_element = CohesiveElement<type>::getNbNodesPerElement();
  UInt nb_nodes_per_sub_element = ElementClass<sub_type>::getNbNodesPerElement();

  UInt nb_points = this->control_points(type, ghost_type).getSize();
  UInt size_of_shapes = CohesiveElement<type>::getShapeSize();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->storage();
  }

  Real * shape_val = shapes.storage();
  Real * u_val     = in_u.storage();
  Vector<Real>::iterator<types::Matrix> uq_it =
    out_uq.begin(nb_points, nb_degree_of_freedom);

  Vector<Real>::const_iterator<types::Matrix> shape =
    shapes.begin(nb_points, nb_nodes_per_element);

  types::Matrix u(nb_nodes_per_element, nb_degree_of_freedom);

  ReduceFunction reduce_function;

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_element;
    if(filter_elements != NULL) {
      shape  = shapes.begin(nb_points, nb_nodes_per_element);
      shape += filter_elem_val[el];
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_sub_element; ++n) {
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	Real u_plus  = u_val[conn_val[el_offset + n] * nb_degree_of_freedom];
	Real u_minus = u_val[conn_val[el_offset + n + nb_nodes_per_sub_element] * nb_degree_of_freedom];
	u(n, d)	= reduce_function(u_plus, u_minus);
      }
    }

    /// @f$ u_q = \textbf{N} * u @f$
    uq_it->mul<false, false>(*shape, u);
    // Math::matrix_matrix(nb_points, nb_degree_of_freedom, nb_nodes_per_element,
    // 			shape, u, uq_val);

    ++uq_it;
    if(filter_elements == NULL) {
      ++shape;
    }
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
