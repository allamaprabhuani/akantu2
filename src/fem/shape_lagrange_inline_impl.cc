/**
 * @file   shape_lagrange_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Feb 15 16:32:44 2011
 *
 * @brief  ShapeLagrange inline implementation
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

__END_AKANTU__

#include "fem.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Vector<Real> & ShapeLagrange<kind>::getShapes(const ElementType & el_type,
							   const GhostType & ghost_type) const {
  return shapes(FEM::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Vector<Real> & ShapeLagrange<kind>::getShapesDerivatives(const ElementType & el_type,
								      const GhostType & ghost_type) const {
  return shapes_derivatives(FEM::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)					\
  setControlPointsByType<type>(control_points, ghost_type);		\
  precomputeShapesOnControlPoints<type>(nodes, ghost_type);		\
  if (ElementClass<type>::getNaturalSpaceDimension() ==			\
      mesh->getSpatialDimension() || kind != _ek_regular)		\
    precomputeShapeDerivativesOnControlPoints<type>(nodes, ghost_type);

template <ElementKind kind>
inline void
ShapeLagrange<kind>::initShapeFunctions(const Vector<Real> & nodes,
					const types::Matrix<Real> & control_points,
					const ElementType & type,
					const GhostType & ghost_type) {
  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
ShapeLagrange<_ek_structural>::initShapeFunctions(__attribute__((unused)) const Vector<Real> & nodes,
						  __attribute__((unused)) const types::Matrix<Real> & control_points,
						  __attribute__((unused)) const ElementType & type,
						  __attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

#undef INIT_SHAPE_FUNCTIONS
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void ShapeLagrange<kind>::
computeShapeDerivativesOnCPointsByElement(const types::Matrix<Real> & node_coords,
					  const types::Matrix<Real> & natural_coords,
					  types::Tensor3<Real> & shapesd) {
  // compute dnds
  types::Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(), natural_coords.cols());
  ElementClass<type>::computeDNDS(natural_coords, dnds);
  // compute dxds
  types::Tensor3<Real> J(node_coords.rows(), natural_coords.rows(), natural_coords.cols());
  ElementClass<type>::computeJMat(dnds, node_coords, J);

  // compute shape derivatives
  ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::inverseMap(const types::RVector & real_coords,
				     UInt elem,
				     types::RVector & natural_coords,
				     const GhostType & ghost_type) const{

  UInt spatial_dimension = mesh->getSpatialDimension();
  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt * elem_val = mesh->getConnectivity(type, ghost_type).values;
  types::RMatrix nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh->extractNodalValuesFromElement(mesh->getNodes(),
				      nodes_coord.storage(),
				      elem_val + elem*nb_nodes_per_element,
				      nb_nodes_per_element,
				      spatial_dimension);

  ElementClass<type>::inverseMap(real_coords,
				 nodes_coord,
				 natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
bool ShapeLagrange<kind>::contains(const types::RVector & real_coords,
			     UInt elem,
			     const GhostType & ghost_type) const{

  UInt spatial_dimension = mesh->getSpatialDimension();
  types::RVector natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  return ElementClass<type>::contains(natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapes(const types::RVector & real_coords,
				  UInt elem,
				  types::RVector & shapes,
				  const GhostType & ghost_type) const{

  UInt spatial_dimension = mesh->getSpatialDimension();
  types::RVector natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  ElementClass<type>::computeShapes(natural_coords, shapes);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
ShapeLagrange<kind>::ShapeLagrange(const Mesh & mesh,
				   const ID & id,
				   const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id),
  shapes("shapes_generic", id),
  shapes_derivatives("shapes_derivatives_generic", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapesOnControlPoints(__attribute__((unused)) const Vector<Real> & nodes,
							  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  types::Matrix<Real> & natural_coords = control_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  UInt nb_element = mesh->getConnectivity(type, ghost_type).getSize();;

  Vector<Real> & shapes_tmp = shapes.alloc(nb_element*nb_points,
					   size_of_shapes,
					   itp_type,
					   ghost_type);

  Vector<Real>::iterator< types::Matrix<Real> > shapes_it =
    shapes_tmp.begin_reinterpret(ElementClass<type>::getNbNodesPerInterpolationElement(), nb_points,
				 nb_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it) {
    types::Matrix<Real> & N = *shapes_it;
    ElementClass<type>::computeShapes(natural_coords,
				      N);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapeDerivativesOnControlPoints(const Vector<Real> & nodes, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  //  Real * coord = mesh->getNodes().values;
  UInt spatial_dimension = mesh->getSpatialDimension();

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();
  UInt size_of_shapesd      = ElementClass<type>::getShapeDerivativesSize();
  types::Matrix<Real> & natural_coords = control_points(type, ghost_type);
  UInt nb_points = natural_coords.cols();

  UInt nb_element = mesh->getConnectivity(type, ghost_type).getSize();
  Vector<Real> & shapes_derivatives_tmp = shapes_derivatives.alloc(nb_element*nb_points,
								   size_of_shapesd,
								   itp_type,
								   ghost_type);

  Vector<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, nodes, x_el,
				  type, ghost_type);

  Real * shapesd_val = shapes_derivatives_tmp.storage();
  Vector<Real>::iterator< types::Matrix<Real> > x_it = x_el.begin(spatial_dimension,
								  nb_nodes_per_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
    types::Matrix<Real> & X = *x_it;
    types::Tensor3<Real> B(shapesd_val,
			   spatial_dimension, nb_nodes_per_element, nb_points);
    computeShapeDerivativesOnCPointsByElement<type>(X,
						    natural_coords,
						    B);

    shapesd_val += size_of_shapesd*nb_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnControlPoints(const Vector<Real> &in_u,
					       Vector<Real> &out_uq,
					       UInt nb_degree_of_freedom,
					       GhostType ghost_type,
					       const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(shapes.exists(itp_type, ghost_type),
		      "No shapes for the type "
		      << shapes.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element = ElementClass<type>::getNbNodesPerInterpolationElement();

  Vector<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, in_u, u_el, type, ghost_type, filter_elements);

  this->interpolateElementalFieldOnControlPoints<type>(u_el, out_uq, ghost_type,
						       shapes(itp_type, ghost_type),
						       filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::gradientOnControlPoints(const Vector<Real> &in_u,
					       Vector<Real> &out_nablauq,
					       UInt nb_degree_of_freedom,
					       GhostType ghost_type,
					       const Vector<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(shapes_derivatives.exists(itp_type, ghost_type),
		      "No shapes derivatives for the type "
		      << shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element  = ElementClass<type>::getNbNodesPerInterpolationElement();

  Vector<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEM::extractNodalToElementField(*mesh, in_u, u_el, type, ghost_type, filter_elements);

  this->gradientElementalFieldOnControlPoints<type>(u_el, out_nablauq, ghost_type,
						    shapes_derivatives(itp_type, ghost_type),
						    filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::fieldTimesShapes(const Vector<Real> & field,
				     Vector<Real> & field_times_shapes,
				     GhostType ghost_type) const {
  field_times_shapes.copy(field);
  field_times_shapes.extendComponentsInterlaced(ElementClass<type>::getShapeSize(), 1);

  UInt nb_element = field_times_shapes.getSize() * ElementClass<type>::getShapeSize();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Real * field_times_shapes_val = field_times_shapes.storage();
  Real * shapes_val = shapes(itp_type, ghost_type).storage();

  /// compute @f$ rho * \varphi_i @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    *field_times_shapes_val++ *= *shapes_val++;
  }
}


/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void ShapeLagrange<kind>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
