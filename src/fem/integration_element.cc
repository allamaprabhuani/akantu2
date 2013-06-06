/**
 * @file   integration_element.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sun Nov 18 00:45:16 2012
 *
 * @brief  Definition of the intagration constants
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
#include "aka_common.hh"
#include "element_class.hh"

using std::sqrt;

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Points                                                                     */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_point, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_positions[]     = {0};
template<> Real GaussIntegrationTypeData<_git_point, 1>::quad_weights[]       = {1.};


/* -------------------------------------------------------------------------- */
/* Segments                                                                   */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_positions[]     = {0.};
template<> Real GaussIntegrationTypeData<_git_segment, 1>::quad_weights[]       = {2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 2>::nb_quadrature_points = 2;
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_positions[]     = {-1./sqrt(3.), 1./sqrt(3.)};
template<> Real GaussIntegrationTypeData<_git_segment, 2>::quad_weights[]       = {1., 1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_segment, 3>::nb_quadrature_points = 3;
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_positions[]     = {-sqrt(3./5.), 0., sqrt(3./5.)};
template<> Real GaussIntegrationTypeData<_git_segment, 3>::quad_weights[]       = {5./9., 8./9., 5./9.};

/* -------------------------------------------------------------------------- */
/* Triangles                                                                  */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_positions[]     = {1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 1>::quad_weights[]       = {1./2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 2>::nb_quadrature_points = 3;
template<> Real GaussIntegrationTypeData<_git_triangle, 2>::quad_positions[]     = {1./6., 1./6.,
										    2./3., 1./6.,
										    1./6., 2./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 2>::quad_weights[]       = {1./6., 1./6., 1./6.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_triangle, 3>::nb_quadrature_points = 4;
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_positions[]     = {1./5., 1./5.,
										    3./5., 1./5.,
										    1./5., 3./5.,
										    1./3., 1./3.};
template<> Real GaussIntegrationTypeData<_git_triangle, 3>::quad_weights[]       = {25./(24.*4.), 25./(24.*4.), 25./(24.*4.), -27/(24.*4.)};

/* -------------------------------------------------------------------------- */
/* Tetrahedrons                                                               */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationTypeData<_git_tetrahedron, 1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_positions[]     = {1./4., 1./4., 1./4.};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 1>::quad_weights[]       = {1./6.};
/* -------------------------------------------------------------------------- */
static const Real tet_2_a = (5. -    std::sqrt(5.))/20.;
static const Real tet_2_b = (5. + 3.*std::sqrt(5.))/20.;
template<> UInt GaussIntegrationTypeData<_git_tetrahedron, 2>::nb_quadrature_points = 4;
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 2>::quad_positions[]     = {tet_2_a, tet_2_a, tet_2_a,
										       tet_2_b, tet_2_a, tet_2_a,
										       tet_2_a, tet_2_b, tet_2_a,
										       tet_2_a, tet_2_a, tet_2_b};
template<> Real GaussIntegrationTypeData<_git_tetrahedron, 2>::quad_weights[]       = {1./24., 1./24., 1./24., 1./24.};



// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationTypeData<_quadrangle_4>::nb_quadrature_points = 4;
// template<> Real GaussIntegrationTypeData<_quadrangle_4>::quad[]               = {-1./sqrt(3.), -1./sqrt(3.),
// 									 1./sqrt(3.), -1./sqrt(3.),
// 									 1./sqrt(3.),  1./sqrt(3.),
// 									 -1./sqrt(3.),  1./sqrt(3.)};
// template<> Real GaussIntegrationTypeData<_quadrangle_4>::weights[]            = {1., 1., 1., 1.};
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationTypeData<_quadrangle_8>::nb_quadrature_points = 9;
// template<> Real GaussIntegrationTypeData<_quadrangle_8>::quad[]               = {0.          ,           0.,
// 									 sqrt(3./5.) ,  sqrt(3./5.),
// 									 -sqrt(3./5.),  sqrt(3./5.),
// 									 -sqrt(3./5.), -sqrt(3./5.),
// 									 sqrt(3./5.) , -sqrt(3./5.),
// 									 0.          ,  sqrt(3./5.),
// 									 -sqrt(3./5.),           0.,
// 									 0.          , -sqrt(3./5.),
// 									 sqrt(3./5.) ,           0.};
// template<> Real GaussIntegrationTypeData<_quadrangle_8    >::weights[] = {64./81.,
// 								  25./81., 25./81., 25./81., 25./81.,
// 								  40./81., 40./81., 40./81., 40./81.};
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationElement<_hexahedron_8>::nb_quadrature_points = 8;
// template<> Real GaussIntegrationElement<_hexahedron_8>::quad[]               = {-1./sqrt(3.), -1./sqrt(3.), -1./sqrt(3.),
// 										1./sqrt(3.) , -1./sqrt(3.), -1./sqrt(3.),
// 										1./sqrt(3.) ,  1./sqrt(3.), -1./sqrt(3.),
// 										-1./sqrt(3.),  1./sqrt(3.), -1./sqrt(3.),
// 										-1./sqrt(3.), -1./sqrt(3.),  1./sqrt(3.),
// 										1./sqrt(3.) , -1./sqrt(3.),  1./sqrt(3.),
// 										1./sqrt(3.) ,  1./sqrt(3.),  1./sqrt(3.),
// 										-1./sqrt(3.),  1./sqrt(3.),  1./sqrt(3.)};
// template<> Real GaussIntegrationElement<_hexahedron_8>::weights[]            = {1., 1., 1., 1.,
// 										1., 1., 1., 1.};


// #if defined(AKANTU_STRUCTURAL_MECHANICS)
// /* -------------------------------------------------------------------------- */
// /* Structural elements                                                        */
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationElement<_bernoulli_beam_2>::nb_quadrature_points = 2;
// template<> Real GaussIntegrationElement<_bernoulli_beam_2>::quad[]               = {-1./sqrt(3.), 0,
// 										    1./sqrt(3.), 0.};
// template<> Real GaussIntegrationElement<_bernoulli_beam_2>::weights[]            = {1., 1.};
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationElement<_bernoulli_beam_3>::nb_quadrature_points  = 3;
// template<> Real GaussIntegrationElement<_bernoulli_beam_3>::quad[]                = {-sqrt(3./5.), 0., 0.,
// 										     0., 0., 0.,
// 										     sqrt(3./5.), 0., 0.};
// template<> Real GaussIntegrationElement<_bernoulli_beam_3>::weights[]             = {5./9., 8./9., 5./9.};
// #endif

// #if defined(AKANTU_COHESIVE_ELEMENT)
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationElement<_cohesive_2d_4>::nb_quadrature_points = 1;
// template<> Real GaussIntegrationElement<_cohesive_2d_4>::quad[]               = {0};
// template<> Real GaussIntegrationElement<_cohesive_2d_4>::weights[]            = {2.};
// /* -------------------------------------------------------------------------- */
// template<> UInt GaussIntegrationElement<_cohesive_2d_6>::nb_quadrature_points = 2;
// template<> Real GaussIntegrationElement<_cohesive_2d_6>::quad[]               = {-1./sqrt(3.), 1./sqrt(3.)};
// template<> Real GaussIntegrationElement<_cohesive_2d_6>::weights[]            = {1., 1.};
// /* -------------------------------------------------------------------------- */
// #endif

__END_AKANTU__
