/**
 * @file   resolution_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Sun Nov 20 2022
 *
 * @brief  Specialization of the resolution class for the penalty method
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

/* -------------------------------------------------------------------------- */
#include "resolution_penalty.hh"
#include "resolution_penalty_tmpl.hh"
#include "penalty_function_linear.hh"
#include "penalty_function_quadratic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

// Instantiate linear penalty as a resolution
using ResolutionPenaltyLinear =
    ResolutionPenalty<Real, PenaltyFunctionLinear<Real>>;
INSTANTIATE_RESOLUTION(penalty_linear, ResolutionPenaltyLinear);

// Instantiate quadratic penalty as a resolution
using ResolutionPenaltyQuadratic =
    ResolutionPenalty<Real, PenaltyFunctionQuadratic<Real>>;
INSTANTIATE_RESOLUTION(penalty_quadratic, ResolutionPenaltyQuadratic);

} // namespace akantu
