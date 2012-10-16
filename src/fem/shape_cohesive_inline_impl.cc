/**
 * @file   shape_cohesive.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
ShapeCohesive<ShapeFunction>::ShapeCohesive(const Mesh & mesh,
					    const ID & id,
					    const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr;
  sstr << id << "sub_shapes_functions";
  sub_type_shape_function = new ShapeFunction(mesh, sstr.str(), memory_id);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
const Vector<Real> & ShapeCohesive<ShapeFunction>::getShapes(const ElementType & type,
							     GhostType ghost_type) const {
  return sub_type_shape_function->getShapes(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
const Vector<Real> & ShapeCohesive<ShapeFunction>::getShapesDerivatives(const ElementType & type,
									GhostType ghost_type) const {
  return sub_type_shape_function->getShapesDerivatives(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
template <ElementType type>
void ShapeCohesive<ShapeFunction>::precomputeShapesOnControlPoints(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);

  Real * natural_coords = this->control_points(type, ghost_type).storage();
  UInt nb_points        = this->control_points(type, ghost_type).getSize();

  UInt size_of_shapes    = CohesiveElement<type>::getShapeSize();

  UInt nb_element = this->mesh->getConnectivity(type, ghost_type).getSize();

  Vector<Real> & shapes_tmp = sub_type_shape_function->shapes.alloc(nb_element*nb_points,
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

  // const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);

  // Real * coord = this->mesh->getNodes().storage();
  //  UInt spatial_dimension = CohesiveElement<type>::getSpatialDimension();

  // UInt nb_nodes_per_element = CohesiveElement<type>::getNbNodesPerElement();
  // UInt nb_nodes_per_sub_element = ElementClass<sub_type>::getNbNodesPerElement();
  UInt size_of_shapesd      = CohesiveElement<type>::getShapeDerivativesSize();
  UInt nb_points            = this->control_points(type, ghost_type).getSize();
  Real * natural_coords     = this->control_points(type, ghost_type).storage();

  // UInt * elem_val = this->mesh->getConnectivity(type, ghost_type).storage();;
  UInt nb_element = this->mesh->getConnectivity(type, ghost_type).getSize();

  Vector<Real> & shapes_derivatives_tmp =
    sub_type_shape_function->shapes_derivatives.alloc(nb_element*nb_points,
						      size_of_shapesd,
						      type,
						      ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.storage();

  for (UInt elem = 0; elem < nb_element; ++elem) {
    CohesiveElement<type>::computeDNDS(natural_coords, nb_points, shapesd_val);
    //    natural_coords += nb_points * spatial_dimension;
    shapesd_val += size_of_shapesd * nb_points;
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

  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);

  AKANTU_DEBUG_ASSERT(sub_type_shape_function->shapes.exists(type, ghost_type),
		      "No shapes for the type "
		      << sub_type_shape_function->shapes.printType(type, ghost_type));

  UInt nb_element = this->mesh->getNbElement(type, ghost_type);
  UInt * conn_val = this->mesh->getConnectivity(type, ghost_type).storage();

  const Vector<Real> & shapes = sub_type_shape_function->shapes(type, ghost_type);

  UInt nb_nodes_per_element = CohesiveElement<type>::getNbNodesPerElement();
  UInt nb_nodes_per_sub_element = ElementClass<sub_type>::getNbNodesPerElement();

  UInt nb_points = this->control_points(type, ghost_type).getSize();
  // UInt size_of_shapes = CohesiveElement<type>::getShapeSize();

  UInt * filter_elem_val = NULL;
  if(filter_elements != NULL) {
    nb_element = filter_elements->getSize();
    filter_elem_val = filter_elements->storage();
  }

  // Real * shape_val = shapes.storage();
  Real * u_val     = in_u.storage();
  Vector<Real>::iterator<types::RMatrix> uq_it =
    out_uq.begin_reinterpret(nb_points, nb_degree_of_freedom,
			     nb_element,
			     nb_points * nb_degree_of_freedom);

  Vector<Real>::const_iterator<types::RMatrix> shape_beginning =
    shapes.begin_reinterpret(nb_points, nb_nodes_per_sub_element,
			     nb_element,
			     nb_points * nb_nodes_per_sub_element);

  types::RMatrix u(nb_nodes_per_sub_element, nb_degree_of_freedom);

  ReduceFunction reduce_function;
  Vector<Real>::const_iterator<types::RMatrix> shape = shape_beginning;

  for (UInt el = 0; el < nb_element; ++el) {
    UInt el_offset = el * nb_nodes_per_sub_element;
    if(filter_elements != NULL) {
      shape = shape_beginning;
      shape += filter_elem_val[el];
      el_offset = filter_elem_val[el] * nb_nodes_per_element;
    }

    for (UInt n = 0; n < nb_nodes_per_sub_element; ++n) {
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	Real u_plus  = u_val[conn_val[el_offset + n] * nb_degree_of_freedom + d];
	Real u_minus = u_val[conn_val[el_offset + n + nb_nodes_per_sub_element] * nb_degree_of_freedom + d];
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

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
template <ElementType type, class ReduceFunction>
void ShapeCohesive<ShapeFunction>::variationOnControlPoints(const Vector<Real> &in_u,
							    Vector<Real> &nablauq,
							    UInt nb_degree_of_freedom,
							    GhostType ghost_type,
							    const Vector<UInt> * filter_elements) const {
  const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);

  AKANTU_DEBUG_ASSERT(sub_type_shape_function->shapes_derivatives.exists(type, ghost_type),
		      "No shapes for the type "
		      << sub_type_shape_function->shapes.printType(type, ghost_type));

  UInt nb_element        = this->mesh->getNbElement(type, ghost_type);
  // UInt spatial_dimension = this->mesh->getSpatialDimension();

  const Vector<UInt> & connectivity = this->mesh->getConnectivity(type, ghost_type);

  const Vector<Real> & shapes_derivatives =
    sub_type_shape_function->shapes_derivatives(type, ghost_type);

  UInt element_dimension        = ElementClass<sub_type>::getSpatialDimension();
  //  UInt nb_nodes_per_element     = CohesiveElement<type>::getNbNodesPerElement();
  UInt nb_nodes_per_sub_element = ElementClass<sub_type>::getNbNodesPerElement();

  UInt nb_points = this->control_points(type, ghost_type).getSize();

  Vector<UInt>::iterator<UInt> filter_elem;
  if(filter_elements != NULL) {
    nb_element  = filter_elements->getSize();
    filter_elem = const_cast<Vector<UInt> *>(filter_elements)->begin();
  }

  Vector<Real>::iterator<types::RMatrix> nablauq_it =
    nablauq.begin(element_dimension, nb_degree_of_freedom);

  Vector<Real>::const_iterator<types::RMatrix> shape_derivative_beginning =
    shapes_derivatives.begin(element_dimension, nb_nodes_per_sub_element);

  types::RMatrix u(nb_nodes_per_sub_element, nb_degree_of_freedom);

  ReduceFunction reduce_function;

  Vector<Real>::const_iterator<types::RMatrix> shape_derivative = shape_derivative_beginning;
  for (UInt el = 0; el < nb_element; ++el) {
    // UInt el_offset = el * nb_nodes_per_element;
    UInt element = el;
    if(filter_elements != NULL) {
      element = *filter_elem;
      shape_derivative = shape_derivative_beginning;
      shape_derivative += *filter_elem;
      ++filter_elem;
    }

    // compute the average/difference of the nodal field loaded from cohesive element
    for (UInt n = 0; n < nb_nodes_per_sub_element; ++n) {
      for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
	Real u_plus  = in_u(connectivity(element, n), d);
	Real u_minus = in_u(connectivity(element, n + nb_nodes_per_sub_element), d);
	u(n, d)	= reduce_function(u_plus, u_minus);
      }
    }

    for (UInt q = 0; q < nb_points; ++q) {
      nablauq_it->mul<false, false>(*shape_derivative, u);
      ++shape_derivative;
      ++nablauq_it;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ShapeFunction>
template <ElementType type, class ReduceFunction>
void ShapeCohesive<ShapeFunction>::computeNormalsOnControlPoints(const Vector<Real> &u,
								 Vector<Real> &normals_u,
								 GhostType ghost_type,
								 const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  // const ElementType sub_type = ElementType(CohesiveElementSubElementType<type>::value);

  UInt nb_element = this->mesh->getNbElement(type, ghost_type);
  UInt nb_points  = this->control_points(type, ghost_type).getSize();
  UInt spatial_dimension = this->mesh->getSpatialDimension();

  if(filter_elements != NULL) {
    nb_element  = filter_elements->getSize();
  }

  Vector<Real> tangents_u(nb_element * nb_points, (spatial_dimension *  (spatial_dimension -1)));

  this->template variationOnControlPoints<type, ReduceFunction>(u,
								tangents_u,
								spatial_dimension,
								ghost_type,
								filter_elements);

  normals_u.resize(nb_points * nb_element);
  Vector<Real>::iterator<types::RVector> normal     = normals_u.begin(spatial_dimension);
  Vector<Real>::iterator<types::RVector> normal_end = normals_u.end(spatial_dimension);

  Real * tangent = tangents_u.storage();

  if(spatial_dimension == 3)
    for (; normal != normal_end; ++normal) {
      Math::vectorProduct3(tangent, tangent+spatial_dimension, normal->storage());

      (*normal) /= normal->norm();

      tangent += spatial_dimension * 2;
    }
  else if(spatial_dimension == 2)
    for (; normal != normal_end; ++normal) {
      types::RVector a1(tangent, spatial_dimension);

      (*normal)(0) = -a1(1);
      (*normal)(1) = a1(0);
      (*normal) /= normal->norm();

      tangent += spatial_dimension;
    }

  AKANTU_DEBUG_OUT();
}

