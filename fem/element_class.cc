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
template<ElementType type> UInt ElementClass<type>::nb_nodes_per_element    = 0;
template<ElementType type> UInt ElementClass<type>::nb_nodes_per_element_p1 = 0;
template<ElementType type> UInt ElementClass<type>::nb_quadrature_points    = 0;
template<ElementType type> UInt ElementClass<type>::spatial_dimension       = 0;
template<ElementType type> ElementType ElementClass<type>::facet_type     = _not_defined;

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_1>::nb_nodes_per_element    = 2;
template<> UInt ElementClass<_line_1>::nb_nodes_per_element_p1 = 2;
template<> UInt ElementClass<_line_1>::nb_quadrature_points    = 1;
template<> UInt ElementClass<_line_1>::spatial_dimension       = 1;
template<> UInt ElementClass<_line_1>::nb_facets               = 2;
template<> ElementType ElementClass<_line_1>::facet_type     = _not_defined;
template<> UInt ElementClass<_line_1>::vec_facet_connectivity[]  = {0,1};
template<> UInt * ElementClass<_line_1>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								  &vec_facet_connectivity[1]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_line_2>::nb_nodes_per_element    = 3;
template<> UInt ElementClass<_line_2>::nb_nodes_per_element_p1 = 2;
template<> UInt ElementClass<_line_2>::nb_quadrature_points    = 2;
template<> UInt ElementClass<_line_2>::spatial_dimension       = 1;
template<> UInt ElementClass<_line_2>::nb_facets               = 2;
template<> ElementType ElementClass<_line_2>::facet_type     = _not_defined;
template<> UInt ElementClass<_line_2>::vec_facet_connectivity[]  = {0,1};
template<> UInt * ElementClass<_line_2>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								  &vec_facet_connectivity[1]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_1>::nb_nodes_per_element    = 3;
template<> UInt ElementClass<_triangle_1>::nb_nodes_per_element_p1 = 3;
template<> UInt ElementClass<_triangle_1>::nb_quadrature_points    = 1;
template<> UInt ElementClass<_triangle_1>::spatial_dimension       = 2;
template<> UInt ElementClass<_triangle_1>::nb_facets               = 3;
template<> ElementType ElementClass<_triangle_1>::facet_type     = _line_1;
template<> UInt ElementClass<_triangle_1>::vec_facet_connectivity[]  = {0,1,1,2,2,0};
template<> UInt * ElementClass<_triangle_1>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								  &vec_facet_connectivity[2],
								  &vec_facet_connectivity[4]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_2>::nb_nodes_per_element    = 6;
template<> UInt ElementClass<_triangle_2>::nb_nodes_per_element_p1 = 3;
template<> UInt ElementClass<_triangle_2>::nb_quadrature_points    = 2;
template<> UInt ElementClass<_triangle_2>::spatial_dimension       = 2;
template<> UInt ElementClass<_triangle_2>::nb_facets               = 3;
template<> ElementType ElementClass<_triangle_2>::facet_type     = _line_2;
template<> UInt ElementClass<_triangle_2>::vec_facet_connectivity[]  = {0,1,3,1,2,4,2,0,5};
template<> UInt * ElementClass<_triangle_2>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								  &vec_facet_connectivity[3],
								  &vec_facet_connectivity[6]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedra_1>::nb_nodes_per_element    = 4;
template<> UInt ElementClass<_tetrahedra_1>::nb_nodes_per_element_p1 = 4;
template<> UInt ElementClass<_tetrahedra_1>::nb_quadrature_points    = 1;
template<> UInt ElementClass<_tetrahedra_1>::spatial_dimension       = 3;
template<> UInt ElementClass<_tetrahedra_1>::nb_facets               = 4;
template<> ElementType ElementClass<_tetrahedra_1>::facet_type     = _triangle_1;
template<> UInt ElementClass<_tetrahedra_1>::vec_facet_connectivity[]  = {0,2,1,1,2,3,2,0,3,0,1,3};
template<> UInt * ElementClass<_tetrahedra_1>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								      &vec_facet_connectivity[3],
								      &vec_facet_connectivity[6],
								      &vec_facet_connectivity[9]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedra_2>::nb_nodes_per_element    = 10;
template<> UInt ElementClass<_tetrahedra_2>::nb_nodes_per_element_p1 = 4;
template<> UInt ElementClass<_tetrahedra_2>::nb_quadrature_points    = 4;
template<> UInt ElementClass<_tetrahedra_2>::spatial_dimension       = 3;
template<> UInt ElementClass<_tetrahedra_2>::nb_facets               = 4;
template<> ElementType ElementClass<_tetrahedra_2>::facet_type     = _triangle_2;
template<> UInt ElementClass<_tetrahedra_2>::vec_facet_connectivity[]  = {0,2,1,6,5,4,1,2,3,5,8,9,2,0,3,8,7,6,0,1,3,4,9,7};
template<> UInt * ElementClass<_tetrahedra_2>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								      &vec_facet_connectivity[3],
								      &vec_facet_connectivity[6],
								      &vec_facet_connectivity[9]};
/* -------------------------------------------------------------------------- */

__END_AKANTU__

