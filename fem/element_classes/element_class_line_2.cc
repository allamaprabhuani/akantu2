/**
 * @file   element_class_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 15 10:28:28 2010
 *
 * @brief  Specialization of the element_class class for the type _line_2
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
template<> ElementClass<_line_2>::ElementClass() {
  nb_nodes_per_element = 3;
  nb_quadrature_points = 2;
  spatial_dimension    = 1;
}

/* -------------------------------------------------------------------------- */


