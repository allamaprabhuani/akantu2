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
/* 2D cohesive elements                                                       */
/* -------------------------------------------------------------------------- */
template<> UInt CohesiveElement<_cohesive_2d_4>::vec_facet_connectivity[]= {0, 1,
									    2, 3};
template<> UInt * CohesiveElement<_cohesive_2d_4>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									    &vec_facet_connectivity[2]};

/* -------------------------------------------------------------------------- */
template<> UInt CohesiveElement<_cohesive_2d_6>::vec_facet_connectivity[]= {0, 1, 2,
									    3, 4, 5};
template<> UInt * CohesiveElement<_cohesive_2d_6>::facet_connectivity[]  = {&vec_facet_connectivity[0],
									    &vec_facet_connectivity[3]};
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/* 3D cohesive elements                                                       */
/* -------------------------------------------------------------------------- */


__END_AKANTU__
