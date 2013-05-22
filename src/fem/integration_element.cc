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
template<> UInt GaussIntegrationElement<_point_1>::nb_quadrature_points = 1;
template<> Real GaussIntegrationElement<_point_1>::quad[]               = {0};
template<> Real GaussIntegrationElement<_point_1>::weights[]            = {1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_segment_2>::nb_quadrature_points = 1;
template<> Real GaussIntegrationElement<_segment_2>::quad[]               = {0};
template<> Real GaussIntegrationElement<_segment_2>::weights[]            = {2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_segment_3>::nb_quadrature_points = 2;
template<> Real GaussIntegrationElement<_segment_3>::quad[]               = {-1./sqrt(3.), 1./sqrt(3.)};
template<> Real GaussIntegrationElement<_segment_3>::weights[]            = {1., 1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_triangle_3>::nb_quadrature_points = 1;
template<> Real GaussIntegrationElement<_triangle_3>::quad[]               = {1./3., 1./3.};
template<> Real GaussIntegrationElement<_triangle_3>::weights[]            = {1./2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_triangle_6>::nb_quadrature_points = 3;
template<> Real GaussIntegrationElement<_triangle_6>::quad[]               = {1./6., 1./6.,
									      2./3., 1./6.,
									      1./6., 2./3.};
template<> Real GaussIntegrationElement<_triangle_6>::weights[]            = {1./6., 1./6., 1./6.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_tetrahedron_4>::nb_quadrature_points = 1;
template<> Real GaussIntegrationElement<_tetrahedron_4>::quad[]               = {1./4., 1./4., 1./4.};
template<> Real GaussIntegrationElement<_tetrahedron_4>::weights[]            = {1./6.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_tetrahedron_10>::nb_quadrature_points = 4;
template<> Real GaussIntegrationElement<_tetrahedron_10>::quad[]               = {0.1381966011250, 0.1381966011250, 0.1381966011250,  // a = (5-sqrt(5))/20
										  0.5854101966250, 0.1381966011250, 0.1381966011250,  // b = (5+3*sqrt(5))/20
										  0.1381966011250, 0.5854101966250, 0.1381966011250,
										  0.1381966011250, 0.1381966011250, 0.5854101966250};
template<> Real GaussIntegrationElement<_tetrahedron_10>::weights[]           = {1./24., 1./24., 1./24., 1./24.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_quadrangle_4>::nb_quadrature_points = 4;
template<> Real GaussIntegrationElement<_quadrangle_4>::quad[]               = {-1./sqrt(3.), -1./sqrt(3.),
										1./sqrt(3.), -1./sqrt(3.),
										1./sqrt(3.),  1./sqrt(3.),
										-1./sqrt(3.),  1./sqrt(3.)};
template<> Real GaussIntegrationElement<_quadrangle_4>::weights[]            = {1., 1., 1., 1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_quadrangle_8>::nb_quadrature_points = 9;
template<> Real GaussIntegrationElement<_quadrangle_8>::quad[]               = {0.          ,           0.,
										sqrt(3./5.) ,  sqrt(3./5.),
										-sqrt(3./5.),  sqrt(3./5.),
										-sqrt(3./5.), -sqrt(3./5.),
										sqrt(3./5.) , -sqrt(3./5.),
										0.          ,  sqrt(3./5.),
										-sqrt(3./5.),           0.,
										0.          , -sqrt(3./5.),
										sqrt(3./5.) ,           0.};
template<> Real GaussIntegrationElement<_quadrangle_8    >::weights[] = {64./81.,
									 25./81., 25./81., 25./81., 25./81.,
									 40./81., 40./81., 40./81., 40./81.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_hexahedron_8>::nb_quadrature_points = 8;
template<> Real GaussIntegrationElement<_hexahedron_8>::quad[]               = {-1./sqrt(3.), -1./sqrt(3.), -1./sqrt(3.),
										1./sqrt(3.) , -1./sqrt(3.), -1./sqrt(3.),
										1./sqrt(3.) ,  1./sqrt(3.), -1./sqrt(3.),
										-1./sqrt(3.),  1./sqrt(3.), -1./sqrt(3.),
										-1./sqrt(3.), -1./sqrt(3.),  1./sqrt(3.),
										1./sqrt(3.) , -1./sqrt(3.),  1./sqrt(3.),
										1./sqrt(3.) ,  1./sqrt(3.),  1./sqrt(3.),
										-1./sqrt(3.),  1./sqrt(3.),  1./sqrt(3.)};
template<> Real GaussIntegrationElement<_hexahedron_8>::weights[]            = {1., 1., 1., 1.,
										1., 1., 1., 1.};


#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
/* Structural elements                                                        */
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_bernoulli_beam_2>::nb_quadrature_points = 2;
template<> Real GaussIntegrationElement<_bernoulli_beam_2>::quad[]               = {-1./sqrt(3.), 0,
										    1./sqrt(3.), 0.};
template<> Real GaussIntegrationElement<_bernoulli_beam_2>::weights[]            = {1., 1.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_bernoulli_beam_3>::nb_quadrature_points  = 3;
template<> Real GaussIntegrationElement<_bernoulli_beam_3>::quad[]                = {-sqrt(3./5.), 0., 0.,
										     0., 0., 0.,
										     sqrt(3./5.), 0., 0.};
template<> Real GaussIntegrationElement<_bernoulli_beam_3>::weights[]             = {5./9., 8./9., 5./9.};
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_cohesive_2d_4>::nb_quadrature_points = 1;
template<> Real GaussIntegrationElement<_cohesive_2d_4>::quad[]               = {0};
template<> Real GaussIntegrationElement<_cohesive_2d_4>::weights[]            = {2.};
/* -------------------------------------------------------------------------- */
template<> UInt GaussIntegrationElement<_cohesive_2d_6>::nb_quadrature_points = 2;
template<> Real GaussIntegrationElement<_cohesive_2d_6>::quad[]               = {-1./sqrt(3.), 1./sqrt(3.)};
template<> Real GaussIntegrationElement<_cohesive_2d_6>::weights[]            = {1., 1.};
/* -------------------------------------------------------------------------- */
#endif

__END_AKANTU__
