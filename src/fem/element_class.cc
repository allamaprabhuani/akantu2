/**
 * @file   element_class.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 20 10:12:44 2010
 *
 * @brief  Common part of element_classes
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<ElementType type> ElementKind ElementClass<type>::kind             = _ek_regular;
template<ElementType type> UInt ElementClass<type>::nb_nodes_per_element    = 0;
template<ElementType type> ElementType ElementClass<type>::p1_element_type  = _not_defined;
template<ElementType type> UInt ElementClass<type>::nb_quadrature_points    = 0;
template<ElementType type> UInt ElementClass<type>::spatial_dimension       = 0;
template<ElementType type> ElementType ElementClass<type>::facet_type       = _not_defined;
template<ElementType type> UInt ElementClass<type>::nb_shape_functions      = 1;


template<> ElementKind ElementClass<_not_defined>::kind             = _ek_regular;
template<> UInt ElementClass<_not_defined>::nb_nodes_per_element    = 0;
template<> ElementType ElementClass<_not_defined>::p1_element_type  = _not_defined;
template<> UInt ElementClass<_not_defined>::nb_quadrature_points    = 0;
template<> UInt ElementClass<_not_defined>::spatial_dimension       = 0;
template<> ElementType ElementClass<_not_defined>::facet_type       = _not_defined;
template<> UInt ElementClass<_not_defined>::nb_shape_functions      = 1;
template<> Real ElementClass<_not_defined>::quad[]                  = {};
template<> UInt ElementClass<_not_defined>::nb_facets               = 0;
template<> UInt * ElementClass<_not_defined>::facet_connectivity[]  = {};
template<> Real ElementClass<_not_defined>::gauss_integration_weights[] = {2.};

/* -------------------------------------------------------------------------- */
/* Element Kind                                                               */
/* -------------------------------------------------------------------------- */
template<> ElementKind ElementClass<_segment_2       >::kind = _ek_regular;
template<> ElementKind ElementClass<_segment_3       >::kind = _ek_regular;
template<> ElementKind ElementClass<_triangle_3      >::kind = _ek_regular;
template<> ElementKind ElementClass<_triangle_6      >::kind = _ek_regular;
template<> ElementKind ElementClass<_tetrahedron_4   >::kind = _ek_regular;
template<> ElementKind ElementClass<_tetrahedron_10  >::kind = _ek_regular;
template<> ElementKind ElementClass<_quadrangle_4    >::kind = _ek_regular;
template<> ElementKind ElementClass<_quadrangle_8    >::kind = _ek_regular;
template<> ElementKind ElementClass<_hexahedron_8    >::kind = _ek_regular;
template<> ElementKind ElementClass<_point           >::kind = _ek_regular;
template<> ElementKind ElementClass<_bernoulli_beam_2>::kind = _ek_regular;
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
template<> ElementKind ElementClass<_cohesive_2d_4   >::kind = _ek_cohesive;
template<> ElementKind ElementClass<_cohesive_2d_6   >::kind = _ek_cohesive;
#endif
/* -------------------------------------------------------------------------- */
/* Regular Elements                                                           */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_point>::nb_shape_functions       = 1;
template<> UInt ElementClass<_point>::nb_nodes_per_element     = 1;
template<> ElementType ElementClass<_point>::p1_element_type   = _point;
template<> UInt ElementClass<_point>::nb_quadrature_points     = 1;
template<> Real ElementClass<_point>::quad[]                   = {0};
template<> ElementType ElementClass<_point>::facet_type        = _not_defined;
template<> UInt ElementClass<_point>::spatial_dimension        = 0;
template<> UInt ElementClass<_point>::nb_facets                = 0;
template<> UInt * ElementClass<_point>::facet_connectivity[]   = {};

/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_segment_2>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_segment_3>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_triangle_3>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_triangle_6>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_tetrahedron_4>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_tetrahedron_10>::nb_shape_functions      = 1;
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
template<> UInt ElementClass<_quadrangle_4>::nb_shape_functions      = 1;
template<> UInt ElementClass<_quadrangle_4>::nb_nodes_per_element    = 4;
template<> ElementType ElementClass<_quadrangle_4>::p1_element_type  = _quadrangle_4;
template<> UInt ElementClass<_quadrangle_4>::nb_quadrature_points    = 4;
//template<> Real ElementClass<_quadrangle_4>::quad[]                  = {0, 0};
template<> Real ElementClass<_quadrangle_4>::quad[]                  = {-1./std::sqrt(3), -1./sqrt(3),
									 1./sqrt(3), -1./sqrt(3),
									 1./sqrt(3),  1./sqrt(3),
									-1./sqrt(3),  1./sqrt(3)};
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
template<> UInt ElementClass<_quadrangle_8>::nb_shape_functions      = 1;
template<> UInt ElementClass<_quadrangle_8>::nb_nodes_per_element    = 8;
template<> ElementType ElementClass<_quadrangle_8>::p1_element_type  = _quadrangle_8;
template<> UInt ElementClass<_quadrangle_8>::nb_quadrature_points    = 9;
template<> Real ElementClass<_quadrangle_8>::quad[]                  = {          0.,           0.,
									 sqrt(3./5.),  sqrt(3./5.),
									-sqrt(3./5.),  sqrt(3./5.),
									-sqrt(3./5.), -sqrt(3./5.),
									 sqrt(3./5.), -sqrt(3./5.),
									          0.,  sqrt(3./5.),
									-sqrt(3./5.),           0.,
										  0., -sqrt(3./5.),
									 sqrt(3./5.),           0.};
template<> UInt ElementClass<_quadrangle_8>::spatial_dimension       = 2;
template<> UInt ElementClass<_quadrangle_8>::nb_facets               = 4;
template<> ElementType ElementClass<_quadrangle_8>::facet_type       = _segment_3;
template<> UInt ElementClass<_quadrangle_8>::vec_facet_connectivity[]= {0, 1, 4,
 									1, 2, 5,
 									2, 3, 6,
									3, 0, 7};
template<> UInt * ElementClass<_quadrangle_8>::facet_connectivity[]  = {vec_facet_connectivity + 0,
 									vec_facet_connectivity + 3,
 									vec_facet_connectivity + 6,
 									vec_facet_connectivity + 9};
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_hexahedron_8>::nb_shape_functions      = 1;
template<> UInt ElementClass<_hexahedron_8>::nb_nodes_per_element    = 8;
template<> ElementType ElementClass<_hexahedron_8>::p1_element_type  = _hexahedron_8;
template<> UInt ElementClass<_hexahedron_8>::nb_quadrature_points    = 8;
template<> Real ElementClass<_hexahedron_8>::quad[]                  = {-1./sqrt(3.), -1./sqrt(3.), -1./sqrt(3.),
 									 1./sqrt(3.), -1./sqrt(3.), -1./sqrt(3.),
									 1./sqrt(3.),  1./sqrt(3.), -1./sqrt(3.),
 									-1./sqrt(3.),  1./sqrt(3.), -1./sqrt(3.),
 									-1./sqrt(3.), -1./sqrt(3.),  1./sqrt(3.),
									 1./sqrt(3.), -1./sqrt(3.),  1./sqrt(3.),
									 1./sqrt(3.),  1./sqrt(3.),  1./sqrt(3.),
 									-1./sqrt(3.),  1./sqrt(3.),  1./sqrt(3.)};
template<> UInt ElementClass<_hexahedron_8>::spatial_dimension       = 3;
template<> UInt ElementClass<_hexahedron_8>::nb_facets               = 6;
template<> ElementType ElementClass<_hexahedron_8>::facet_type       = _quadrangle_4;
template<> UInt ElementClass<_hexahedron_8>::vec_facet_connectivity[]= {0, 1, 2, 3,
									0, 1, 5, 4,
									1, 2, 6, 5,
									2, 3, 7, 6,
									3, 0, 4, 7,
									4, 5, 6, 7};
template<> UInt * ElementClass<_hexahedron_8>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									&vec_facet_connectivity[4],
									&vec_facet_connectivity[8],
									&vec_facet_connectivity[12],
									&vec_facet_connectivity[16],
									&vec_facet_connectivity[20]};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Structural elements                                                        */
/* -------------------------------------------------------------------------- */
template<> UInt ElementClass<_bernoulli_beam_2>::nb_nodes_per_element    = 2;
template<> ElementType ElementClass<_bernoulli_beam_2>::p1_element_type  = _bernoulli_beam_2;
template<> UInt ElementClass<_bernoulli_beam_2>::nb_quadrature_points    = 2;
template<> Real ElementClass<_bernoulli_beam_2>::quad[]                  = {-1/sqrt(3), 0,
									     1/sqrt(3), 0};
template<> UInt ElementClass<_bernoulli_beam_2>::spatial_dimension       = 2;
template<> UInt ElementClass<_bernoulli_beam_2>::nb_facets               = 2;
template<> ElementType ElementClass<_bernoulli_beam_2>::facet_type       = _point;
template<> UInt ElementClass<_bernoulli_beam_2>::vec_facet_connectivity[]= {0,
									    1};
template<> UInt * ElementClass<_bernoulli_beam_2>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									    &vec_facet_connectivity[1]};
template<> UInt ElementClass<_bernoulli_beam_2>::nb_shape_functions      = 5;
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Gauss integration                                                          */
/* -------------------------------------------------------------------------- */
template<> Real ElementClass<_segment_2       >::gauss_integration_weights[] = {2.};
template<> Real ElementClass<_segment_3       >::gauss_integration_weights[] = {1., 1.};
template<> Real ElementClass<_triangle_3      >::gauss_integration_weights[] = {1./2.};
template<> Real ElementClass<_triangle_6      >::gauss_integration_weights[] = {1./6., 1./6., 1./6.};
template<> Real ElementClass<_tetrahedron_4   >::gauss_integration_weights[] = {1./6.};
template<> Real ElementClass<_tetrahedron_10  >::gauss_integration_weights[] = {1./24., 1./24., 1./24., 1./24.};
template<> Real ElementClass<_quadrangle_4    >::gauss_integration_weights[] = {1., 1., 1., 1.};
template<> Real ElementClass<_quadrangle_8    >::gauss_integration_weights[] = {64./81.,
										25./81., 25./81., 25./81., 25./81.,
										40./81., 40./81., 40./81., 40./81.};
template<> Real ElementClass<_hexahedron_8    >::gauss_integration_weights[] = {1., 1., 1., 1.,
					      				      1., 1., 1., 1.};
template<> Real ElementClass<_point           >::gauss_integration_weights[] = {1.};
template<> Real ElementClass<_bernoulli_beam_2>::gauss_integration_weights[] = {1., 1.};


__END_AKANTU__

