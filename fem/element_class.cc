/**
 * @file   element_class.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 20 10:12:44 2010
 *
 * @brief  Common part of element_classes
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<ElementType type> UInt ElementClass<type>::nb_nodes_per_element    = 0;
template<ElementType type> ElementType ElementClass<type>::p1_element_type  = _not_defined;
template<ElementType type> UInt ElementClass<type>::nb_quadrature_points    = 0;
template<ElementType type> UInt ElementClass<type>::spatial_dimension       = 0;
template<ElementType type> ElementType ElementClass<type>::facet_type       = _not_defined;

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_point>::nb_nodes_per_element     = 1;
template<> ElementType ElementClass<_point>::p1_element_type   = _point;
template<> ElementType ElementClass<_point>::facet_type        = _not_defined;
template<> UInt ElementClass<_point>::spatial_dimension        = 0;
template<> UInt ElementClass<_point>::nb_facets                = 0;
template<> UInt * ElementClass<_point>::facet_connectivity[]   = {};

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_segment_2>::nb_nodes_per_element    = 2;
template<> ElementType ElementClass<_segment_2>::p1_element_type  = _segment_2;
template<> UInt ElementClass<_segment_2>::nb_quadrature_points    = 1;
template<> Real ElementClass<_segment_2>::quad[]                  = {0};
template<> UInt ElementClass<_segment_2>::spatial_dimension       = 1;
template<> UInt ElementClass<_segment_2>::nb_facets               = 2;
template<> ElementType ElementClass<_segment_2>::facet_type       = _point;
template<> UInt ElementClass<_segment_2>::vec_facet_connectivity[]= {0,
								     1};
template<> UInt * ElementClass<_segment_2>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								     &vec_facet_connectivity[1]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_segment_3>::nb_nodes_per_element    = 3;
template<> ElementType ElementClass<_segment_3>::p1_element_type  = _segment_2;
template<> UInt ElementClass<_segment_3>::nb_quadrature_points    = 2;
template<> Real ElementClass<_segment_3>::quad[]                  = {-1./sqrt(3.), 1./sqrt(3.)};
template<> UInt ElementClass<_segment_3>::spatial_dimension       = 1;
template<> UInt ElementClass<_segment_3>::nb_facets               = 2;
template<> ElementType ElementClass<_segment_3>::facet_type       = _point;
template<> UInt ElementClass<_segment_3>::vec_facet_connectivity[]= {0,
								     1};
template<> UInt * ElementClass<_segment_3>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								     &vec_facet_connectivity[1]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_3>::nb_nodes_per_element    = 3;
template<> ElementType ElementClass<_triangle_3>::p1_element_type  = _triangle_3;
template<> UInt ElementClass<_triangle_3>::nb_quadrature_points    = 1;
template<> Real ElementClass<_triangle_3>::quad[]                  = {1./3., 1./3.};
template<> UInt ElementClass<_triangle_3>::spatial_dimension       = 2;
template<> UInt ElementClass<_triangle_3>::nb_facets               = 3;
template<> ElementType ElementClass<_triangle_3>::facet_type       = _segment_2;
template<> UInt ElementClass<_triangle_3>::vec_facet_connectivity[]= {0, 1,
								      1, 2,
								      2, 0};
template<> UInt * ElementClass<_triangle_3>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								      &vec_facet_connectivity[2],
								      &vec_facet_connectivity[4]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_triangle_6>::nb_nodes_per_element    = 6;
template<> ElementType ElementClass<_triangle_6>::p1_element_type  = _triangle_3;
template<> UInt ElementClass<_triangle_6>::nb_quadrature_points    = 3;
template<> Real ElementClass<_triangle_6>::quad[]                  = {1./6., 1./6.,
								      2./3., 1./6.,
								      1./6., 2./3.};
template<> UInt ElementClass<_triangle_6>::spatial_dimension       = 2;
template<> UInt ElementClass<_triangle_6>::nb_facets               = 3;
template<> ElementType ElementClass<_triangle_6>::facet_type       = _segment_3;
template<> UInt ElementClass<_triangle_6>::vec_facet_connectivity[]= {0, 1, 3,
								      1, 2, 4,
								      2, 0, 5};
template<> UInt * ElementClass<_triangle_6>::facet_connectivity[]  = {&vec_facet_connectivity[0],
								      &vec_facet_connectivity[3],
								      &vec_facet_connectivity[6]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedron_4>::nb_nodes_per_element    = 4;
template<> ElementType ElementClass<_tetrahedron_4>::p1_element_type  = _tetrahedron_4;
template<> UInt ElementClass<_tetrahedron_4>::nb_quadrature_points    = 1;
template<> Real ElementClass<_tetrahedron_4>::quad[]                  = {1./4., 1./4., 1./4.};
template<> UInt ElementClass<_tetrahedron_4>::spatial_dimension       = 3;
template<> UInt ElementClass<_tetrahedron_4>::nb_facets               = 4;
template<> ElementType ElementClass<_tetrahedron_4>::facet_type       = _triangle_3;
template<> UInt ElementClass<_tetrahedron_4>::vec_facet_connectivity[]= {0, 2, 1,
									 1, 2, 3,
									 2, 0, 3,
									 0, 1, 3};
template<> UInt * ElementClass<_tetrahedron_4>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									 &vec_facet_connectivity[3],
									 &vec_facet_connectivity[6],
									 &vec_facet_connectivity[9]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_tetrahedron_10>::nb_nodes_per_element    = 10;
template<> ElementType ElementClass<_tetrahedron_10>::p1_element_type  = _tetrahedron_4;
template<> UInt ElementClass<_tetrahedron_10>::nb_quadrature_points    = 4;
template<> Real ElementClass<_tetrahedron_10>::quad[] = {0.1381966011250, 0.1381966011250, 0.1381966011250,  // a = (5-sqrt(5))/20
                                                         0.5854101966250, 0.1381966011250, 0.1381966011250,  // b = (5+3*sqrt(5))/20
                                                         0.1381966011250, 0.5854101966250, 0.1381966011250,
                                                         0.1381966011250, 0.1381966011250, 0.5854101966250};
template<> UInt ElementClass<_tetrahedron_10>::spatial_dimension       = 3;
template<> UInt ElementClass<_tetrahedron_10>::nb_facets               = 4;
template<> ElementType ElementClass<_tetrahedron_10>::facet_type       = _triangle_6;
template<> UInt ElementClass<_tetrahedron_10>::vec_facet_connectivity[]= {0, 2, 1, 6, 5, 4,
									  1, 2, 3, 5, 9, 8,
									  2, 0, 3, 6, 7, 9,
									  0, 1, 3, 4, 8, 7};
template<> UInt * ElementClass<_tetrahedron_10>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									  &vec_facet_connectivity[6],
									  &vec_facet_connectivity[12],
									  &vec_facet_connectivity[18]};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_quadrangle_4>::nb_nodes_per_element    = 4;
template<> ElementType ElementClass<_quadrangle_4>::p1_element_type  = _quadrangle_4;
template<> UInt ElementClass<_quadrangle_4>::nb_quadrature_points    = 1;
template<> Real ElementClass<_quadrangle_4>::quad[]                  = {0, 0};
template<> UInt ElementClass<_quadrangle_4>::spatial_dimension       = 2;
template<> UInt ElementClass<_quadrangle_4>::nb_facets               = 4;
template<> ElementType ElementClass<_quadrangle_4>::facet_type       = _segment_2;
template<> UInt ElementClass<_quadrangle_4>::vec_facet_connectivity[]= {0, 1,
									1, 2,
									2, 3,
                                                                        3, 0};
template<> UInt * ElementClass<_quadrangle_4>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									&vec_facet_connectivity[2],
									&vec_facet_connectivity[4],
									&vec_facet_connectivity[6]};
/* -------------------------------------------------------------------------- */

__END_AKANTU__

