/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _tetrahedra_2
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<> ElementClass<_tetrahedra_2>::ElementClass() {
  nb_nodes_per_element = 10;
  nb_quadrature_points = 4;
  spatial_dimension    = 3;
}

/* -------------------------------------------------------------------------- */


