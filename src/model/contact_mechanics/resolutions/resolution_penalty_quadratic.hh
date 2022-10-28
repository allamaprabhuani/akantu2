/**
 * @file   resolution_penalty_quadratic.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Oct 21 2022
 *
 * @brief  Quadratic Penalty Resolution for Contact Mechanics Model
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "resolution_penalty.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH__
#define __AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH__

namespace akantu {

// Define function and derivative for quadratic penalty method
template <typename T> T quadraticPenetrationFunction(T var) {
  return var*var + var;
}

template <typename T> T quadraticPenetrationDerivative(T var) {
  return 2.0*var + 1.0;
}

// Instantiate quadratic penalty as a resolution
INSTANTIATE_RESOLUTION(penalty_quadratic, ResolutionPenaltyTemplate, Real,
                       quadraticPenetrationFunction, quadraticPenetrationDerivative);

} // namespace akantu

#endif /* __AKANTU_RESOLUTION_PENALTY_QUADRATIC_HH__ */
