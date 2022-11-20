/**
 * @file   penalty_function_quadratic.hh
 *
 * @author Fabio Matti <fabio.matti@epfl.ch>
 *
 * @date creation: Fri Oct 28 2022
 * @date last modification: Sun Nov 20 2022
 *
 * @brief  Quadratic penalty function
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_PENALTY_FUNCTION_QUADRATIC_HH__
#define __AKANTU_PENALTY_FUNCTION_QUADRATIC_HH__

/* -------------------------------------------------------------------------- */

template <typename T> class PenaltyFunctionQuadratic {
public:
  T operator()(const T &var) const { return var * var + var; }
  T derivative(const T &var) const { return 2.0 * var + 1.0; }
};

#endif /* __AKANTU_PENALTY_FUNCTION_QUADRATIC_HH__  */
