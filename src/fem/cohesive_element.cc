/**
 * @file   cohesive_element.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Feb 03 16:28:03 2012
 *
 * @brief  CohesiveElement implementation
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
#include "cohesive_element.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_cohesive_2d_4>::spatial_dimension    = 2;
template<> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_nodes_per_element = 4;
template<> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_facets            = 2;
template<> UInt GeometricalElement<_gt_cohesive_2d_4>::nb_nodes_per_facet   = 2;
template<> UInt GeometricalElement<_gt_cohesive_2d_4>::facet_connectivity[] = {0, 2,
									       1, 3};

template<> ElementType ElementClass<_cohesive_2d_4>::facet_type = _segment_2;
template<> ElementType ElementClass<_cohesive_2d_4>::p1_type    = _cohesive_2d_4;

/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_cohesive_2d_6>::spatial_dimension    = 2;
template<> UInt GeometricalElement<_gt_cohesive_2d_6>::nb_nodes_per_element = 6;
template<> UInt GeometricalElement<_gt_cohesive_2d_6>::nb_facets            = 2;
template<> UInt GeometricalElement<_gt_cohesive_2d_6>::nb_nodes_per_facet   = 3;
template<> UInt GeometricalElement<_gt_cohesive_2d_6>::facet_connectivity[] = {0, 3,
									       1, 4,
									       2, 5};

template<> ElementType ElementClass<_cohesive_2d_6>::facet_type = _segment_3;
template<> ElementType ElementClass<_cohesive_2d_6>::p1_type    = _cohesive_2d_4;
/* -------------------------------------------------------------------------- */
template<> UInt GeometricalElement<_gt_cohesive_1d_2>::spatial_dimension    = 1;
template<> UInt GeometricalElement<_gt_cohesive_1d_2>::nb_nodes_per_element = 2;
template<> UInt GeometricalElement<_gt_cohesive_1d_2>::nb_facets            = 2;
template<> UInt GeometricalElement<_gt_cohesive_1d_2>::nb_nodes_per_facet   = 1;
template<> UInt GeometricalElement<_gt_cohesive_1d_2>::facet_connectivity[] = {0,
									       1};

template<> ElementType ElementClass<_cohesive_1d_2>::facet_type = _point_1;
template<> ElementType ElementClass<_cohesive_1d_2>::p1_type    = _cohesive_1d_2;
/* -------------------------------------------------------------------------- */
__END_AKANTU__
