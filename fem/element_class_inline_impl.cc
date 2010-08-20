/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Implementation of the inline functions of the class element_class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<ElementType type> inline UInt ElementClass<type>::getShapeSize() {
  return nb_quadrature_points * nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
template<ElementType type> inline UInt ElementClass<type>::getShapeDerivativesSize() {
  return nb_quadrature_points * nb_nodes_per_element * spatial_dimension;
}

/* -------------------------------------------------------------------------- */
template<ElementType type> inline UInt ElementClass<type>::getJacobianSize() {
  return nb_quadrature_points;
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
void ElementClass<type>::shapeFunctions(__attribute__ ((unused)) const Real * coord,
					__attribute__ ((unused)) Real * shape,
					__attribute__ ((unused)) Real * dshape,
					__attribute__ ((unused)) Real * jacobians) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline Real ElementClass<type>::getInradius(__attribute__ ((unused)) const Real * coord) {
  AKANTU_DEBUG_ERROR("Function not implemented for type : " << type);
}

/* -------------------------------------------------------------------------- */

#include "element_classes/element_class_line_1.cc"
// #include "element_classes/element_class_line_2.cc"
#include "element_classes/element_class_triangle_1.cc"
// #include "element_classes/element_class_triangle_2.cc"
#include "element_classes/element_class_tetrahedra_1.cc"
// #include "element_classes/element_class_tetrahedra_2.cc"
