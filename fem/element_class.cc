/**
 * @file   element_class.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 20 10:12:44 2010
 *
 * @brief  Common part of element_classes
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<ElementType type> UInt ElementClass<type>::nb_nodes_per_element = 0;
template<ElementType type> UInt ElementClass<type>::nb_quadrature_points = 0;
template<ElementType type> UInt ElementClass<type>::spatial_dimension    = 0;

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_1>::nb_nodes_per_element = 2;
template<> UInt ElementClass<_line_1>::nb_quadrature_points = 1;
template<> UInt ElementClass<_line_1>::spatial_dimension    = 1;
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_2>::nb_nodes_per_element = 3;
template<> UInt ElementClass<_line_2>::nb_quadrature_points = 2;
template<> UInt ElementClass<_line_2>::spatial_dimension    = 1;
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_1>::nb_nodes_per_element = 3;
template<> UInt ElementClass<_triangle_1>::nb_quadrature_points = 1;
template<> UInt ElementClass<_triangle_1>::spatial_dimension    = 2;
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_2>::nb_nodes_per_element = 6;
template<> UInt ElementClass<_triangle_2>::nb_quadrature_points = 2;
template<> UInt ElementClass<_triangle_2>::spatial_dimension    = 2;
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedra_1>::nb_nodes_per_element = 4;
template<> UInt ElementClass<_tetrahedra_1>::nb_quadrature_points = 1;
template<> UInt ElementClass<_tetrahedra_1>::spatial_dimension    = 3;
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedra_2>::nb_nodes_per_element = 10;
template<> UInt ElementClass<_tetrahedra_2>::nb_quadrature_points = 4;
template<> UInt ElementClass<_tetrahedra_2>::spatial_dimension    = 3;
/* -------------------------------------------------------------------------- */

__END_AKANTU__

