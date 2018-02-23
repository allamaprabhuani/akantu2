/**
 * @file   non_linear_solver_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Feb 16 10:37:01 2016
 *
 * @brief  Include for the default non linear solvers
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
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__
#define __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__

#if defined(AKANTU_IMPLICIT)
#include "non_linear_solver_linear.hh"
#include "non_linear_solver_newton_raphson.hh"
#endif

#include "non_linear_solver_lumped.hh"

#endif /* __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__ */
