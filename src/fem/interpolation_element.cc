/**
 * @file   interpolation_element.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Jul 20 23:40:43 2010
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
template<> UInt InterpolationElement<_itp_not_defined>::nb_nodes_per_element    = 0;
template<> UInt InterpolationElement<_itp_not_defined>::natural_space_dimension = 0;

/* -------------------------------------------------------------------------- */
/* Regular Elements                                                           */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_point_1>::nb_nodes_per_element    = 1;
template<> UInt InterpolationElement<_itp_lagrange_point_1>::natural_space_dimension = 0;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_segment_2>::nb_nodes_per_element    = 2;
template<> UInt InterpolationElement<_itp_lagrange_segment_2>::natural_space_dimension = 1;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_segment_3>::nb_nodes_per_element    = 3;
template<> UInt InterpolationElement<_itp_lagrange_segment_3>::natural_space_dimension = 1;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_triangle_3>::nb_nodes_per_element    = 3;
template<> UInt InterpolationElement<_itp_lagrange_triangle_3>::natural_space_dimension = 2;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_triangle_6>::nb_nodes_per_element    = 6;
template<> UInt InterpolationElement<_itp_lagrange_triangle_6>::natural_space_dimension = 2;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_tetrahedron_4>::nb_nodes_per_element    = 4;
template<> UInt InterpolationElement<_itp_lagrange_tetrahedron_4>::natural_space_dimension = 3;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_tetrahedron_10>::nb_nodes_per_element    = 10;
template<> UInt InterpolationElement<_itp_lagrange_tetrahedron_10>::natural_space_dimension = 3;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_quadrangle_4>::nb_nodes_per_element    = 4;
template<> UInt InterpolationElement<_itp_lagrange_quadrangle_4>::natural_space_dimension = 2;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_serendip_quadrangle_8>::nb_nodes_per_element    = 8;
template<> UInt InterpolationElement<_itp_serendip_quadrangle_8>::natural_space_dimension = 2;
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_lagrange_hexahedron_8>::nb_nodes_per_element  = 8;
template<> UInt InterpolationElement<_itp_lagrange_hexahedron_8>::natural_space_dimension = 3;
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_STRUCTURAL_MECHANICS)
/* -------------------------------------------------------------------------- */
/* Structural elements                                                        */
/* -------------------------------------------------------------------------- */
template<> UInt InterpolationElement<_itp_bernoulli_beam>::nb_nodes_per_element    = 2;
template<> UInt InterpolationElement<_itp_bernoulli_beam>::natural_space_dimension = 1;
template<> UInt InterpolationElement<_itp_bernoulli_beam>::nb_shape_functions = 5;
/* -------------------------------------------------------------------------- */
#endif

__END_AKANTU__
