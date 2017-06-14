/**
 * @file   shape_lagrange_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 27 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  ShapeLagrange inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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

#include "fe_engine.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Array<Real> &
ShapeLagrange<kind>::getShapes(const ElementType & el_type,
                               const GhostType & ghost_type) const {
  return shapes(FEEngine::getInterpolationType(el_type), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline const Array<Real> &
ShapeLagrange<kind>::getShapesDerivatives(const ElementType & el_type,
                                          const GhostType & ghost_type) const {
  return shapes_derivatives(FEEngine::getInterpolationType(el_type),
                            ghost_type);
}

/* -------------------------------------------------------------------------- */
#define INIT_SHAPE_FUNCTIONS(type)                                             \
  setIntegrationPointsByType<type>(integration_points, ghost_type);            \
  precomputeShapesOnIntegrationPoints<type>(nodes, ghost_type);                \
  if (ElementClass<type>::getNaturalSpaceDimension() ==                        \
          mesh.getSpatialDimension() ||                                        \
      kind != _ek_regular)                                                     \
    precomputeShapeDerivativesOnIntegrationPoints<type>(nodes, ghost_type);

template <ElementKind kind>
inline void ShapeLagrange<kind>::initShapeFunctions(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    const ElementType & type, const GhostType & ghost_type) {
  AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(INIT_SHAPE_FUNCTIONS);
}

#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
template <>
inline void ShapeLagrange<_ek_structural>::initShapeFunctions(
    __attribute__((unused)) const Array<Real> & nodes,
    __attribute__((unused)) const Matrix<Real> & integration_points,
    __attribute__((unused)) const ElementType & type,
    __attribute__((unused)) const GhostType & ghost_type) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}
#endif

#undef INIT_SHAPE_FUNCTIONS
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void ShapeLagrange<kind>::computeShapeDerivativesOnCPointsByElement(
    const Matrix<Real> & node_coords, const Matrix<Real> & natural_coords,
    Tensor3<Real> & shapesd) const {
  AKANTU_DEBUG_IN();

  // compute dnds
  Tensor3<Real> dnds(node_coords.rows(), node_coords.cols(),
                     natural_coords.cols());
  ElementClass<type>::computeDNDS(natural_coords, dnds);
  // compute dxds
  Tensor3<Real> J(node_coords.rows(), natural_coords.rows(),
                  natural_coords.cols());
  ElementClass<type>::computeJMat(dnds, node_coords, J);

  // compute shape derivatives
  ElementClass<type>::computeShapeDerivatives(J, dnds, shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::inverseMap(const Vector<Real> & real_coords,
                                     UInt elem, Vector<Real> & natural_coords,
                                     const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(mesh.getNodes(), nodes_coord.storage(),
                                     elem_val + elem * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  ElementClass<type>::inverseMap(real_coords, nodes_coord, natural_coords);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
bool ShapeLagrange<kind>::contains(const Vector<Real> & real_coords, UInt elem,
                                   const GhostType & ghost_type) const {

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  return ElementClass<type>::contains(natural_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolate(const Vector<Real> & real_coords,
                                      UInt elem,
                                      const Matrix<Real> & nodal_values,
                                      Vector<Real> & interpolated,
                                      const GhostType & ghost_type) const {
  UInt nb_shapes = ElementClass<type>::getShapeSize();
  Vector<Real> shapes(nb_shapes);
  computeShapes<type>(real_coords, elem, shapes, ghost_type);
  ElementClass<type>::interpolate(nodal_values, shapes, interpolated);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapes(const Vector<Real> & real_coords,
                                        UInt elem, Vector<Real> & shapes,
                                        const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> natural_coords(spatial_dimension);

  inverseMap<type>(real_coords, elem, natural_coords, ghost_type);
  ElementClass<type>::computeShapes(natural_coords, shapes);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapeDerivatives(
    const Matrix<Real> & real_coords, UInt elem, Tensor3<Real> & shapesd,
    const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_points = real_coords.cols();
  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == shapesd.size(0) &&
                          nb_nodes_per_element == shapesd.size(1),
                      "Shape size doesn't match");
  AKANTU_DEBUG_ASSERT(nb_points == shapesd.size(2),
                      "Number of points doesn't match shapes size");

  Matrix<Real> natural_coords(spatial_dimension, nb_points);

  // Creates the matrix of natural coordinates
  for (UInt i = 0; i < nb_points; i++) {
    Vector<Real> real_point = real_coords(i);
    Vector<Real> natural_point = natural_coords(i);

    inverseMap<type>(real_point, elem, natural_point, ghost_type);
  }

  UInt * elem_val = mesh.getConnectivity(type, ghost_type).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(mesh.getNodes(), nodes_coord.storage(),
                                     elem_val + elem * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

  computeShapeDerivativesOnCPointsByElement<type>(nodes_coord, natural_coords,
                                                  shapesd);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
ShapeLagrange<kind>::ShapeLagrange(const Mesh & mesh, const ID & id,
                                   const MemoryID & memory_id)
    : ShapeFunctions(mesh, id, memory_id),
      shapes("shapes_generic", id, memory_id),
      shapes_derivatives("shapes_derivatives_generic", id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapesOnIntegrationPoints(
    __attribute__((unused)) const Array<Real> & nodes,
    const Matrix<Real> & integration_points, Array<Real> & shapes,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();

  shapes.resize(nb_element * nb_points);

#if !defined(AKANTU_NDEBUG)
  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  AKANTU_DEBUG_ASSERT(shapes.getNbComponent() == size_of_shapes,
                      "The shapes array does not have the correct "
                          << "number of component");
#endif

  Array<Real>::matrix_iterator shapes_it = shapes.begin_reinterpret(
      ElementClass<type>::getNbNodesPerInterpolationElement(), nb_points,
      nb_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++shapes_it) {
    Matrix<Real> & N = *shapes_it;
    ElementClass<type>::computeShapes(integration_points, N);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  Matrix<Real> & natural_coords = integration_points(type, ghost_type);
  UInt size_of_shapes = ElementClass<type>::getShapeSize();

  Array<Real> & shapes_tmp =
      shapes.alloc(0,
                   size_of_shapes, itp_type, ghost_type);

  this->computeShapesOnIntegrationPoints<type>(nodes, natural_coords,
                                               shapes_tmp, ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::computeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const Matrix<Real> & integration_points,
    Array<Real> & shape_derivatives, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  UInt nb_points = integration_points.cols();
  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();

  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();
  AKANTU_DEBUG_ASSERT(shape_derivatives.getNbComponent() == size_of_shapesd,
                      "The shapes_derivatives array does not have the correct "
                          << "number of component");
  shape_derivatives.resize(nb_element * nb_points);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Real * shapesd_val = shape_derivatives.storage();
  Array<Real>::matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
    Matrix<Real> & X = *x_it;
    Tensor3<Real> B(shapesd_val, spatial_dimension, nb_nodes_per_element,
                    nb_points);
    computeShapeDerivativesOnCPointsByElement<type>(X, integration_points, B);

    shapesd_val += size_of_shapesd * nb_points;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::precomputeShapeDerivativesOnIntegrationPoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Matrix<Real> & natural_coords = integration_points(type, ghost_type);

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;

  UInt size_of_shapesd = ElementClass<type>::getShapeDerivativesSize();

  Array<Real> & shapes_derivatives_tmp =
      shapes_derivatives.alloc(0,
                               size_of_shapesd, itp_type, ghost_type);

  this->computeShapeDerivativesOnIntegrationPoints<type>(
      nodes, natural_coords, shapes_derivatives_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    const Array<Real> & shapes, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->interpolateElementalFieldOnIntegrationPoints<type>(
      u_el, out_uq, ghost_type, shapes, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::interpolateOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_uq, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(shapes.exists(itp_type, ghost_type),
                      "No shapes for the type "
                          << shapes.printType(itp_type, ghost_type));

  this->interpolateOnIntegrationPoints<type>(in_u, out_uq, nb_degree_of_freedom,
                                             shapes(itp_type, ghost_type),
                                             ghost_type, filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::gradientOnIntegrationPoints(
    const Array<Real> & in_u, Array<Real> & out_nablauq,
    UInt nb_degree_of_freedom, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  AKANTU_DEBUG_ASSERT(
      shapes_derivatives.exists(itp_type, ghost_type),
      "No shapes derivatives for the type "
          << shapes_derivatives.printType(itp_type, ghost_type));

  UInt nb_nodes_per_element =
      ElementClass<type>::getNbNodesPerInterpolationElement();

  Array<Real> u_el(0, nb_degree_of_freedom * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, in_u, u_el, type, ghost_type,
                                       filter_elements);

  this->gradientElementalFieldOnIntegrationPoints<type>(
      u_el, out_nablauq, ghost_type, shapes_derivatives(itp_type, ghost_type),
      filter_elements);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void ShapeLagrange<kind>::fieldTimesShapes(const Array<Real> & field,
                                           Array<Real> & field_times_shapes,
                                           GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  field_times_shapes.resize(field.getSize());

  UInt size_of_shapes = ElementClass<type>::getShapeSize();
  InterpolationType itp_type = ElementClassProperty<type>::interpolation_type;
  UInt nb_degree_of_freedom = field.getNbComponent();

  const Array<Real> & shape = shapes(itp_type, ghost_type);

  Array<Real>::const_matrix_iterator field_it =
      field.begin(nb_degree_of_freedom, 1);
  Array<Real>::const_matrix_iterator shapes_it = shape.begin(1, size_of_shapes);

  Array<Real>::matrix_iterator it =
      field_times_shapes.begin(nb_degree_of_freedom, size_of_shapes);
  Array<Real>::matrix_iterator end =
      field_times_shapes.end(nb_degree_of_freedom, size_of_shapes);

  for (; it != end; ++it, ++field_it, ++shapes_it) {
    it->mul<false, false>(*field_it, *shapes_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
void ShapeLagrange<kind>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "Shapes Lagrange [" << std::endl;
  ShapeFunctions::printself(stream, indent + 1);
  shapes.printself(stream, indent + 1);
  shapes_derivatives.printself(stream, indent + 1);
  stream << space << "]" << std::endl;
}
