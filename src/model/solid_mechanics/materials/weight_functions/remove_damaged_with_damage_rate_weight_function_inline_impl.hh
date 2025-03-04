/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "remove_damaged_with_damage_rate_weight_function.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void RemoveDamagedWithDamageRateWeightFunction::init() {
  this->damage_with_damage_rate =
      &(this->manager.registerWeightFunctionInternal("damage-rate"));
}

/* -------------------------------------------------------------------------- */
inline Real RemoveDamagedWithDamageRateWeightFunction::operator()(
    Real r, const IntegrationPoint & q1, const IntegrationPoint & q2) {

  /// compute the weight
  auto quad = q2.global_num;

  if (q1.global_num == quad) {
    return 1.;
  }

  Array<Real> & dam_array =
      (*this->damage_with_damage_rate)(q2.type, q2.ghost_type);
  Real D = dam_array(quad);
  Real w = 0.;
  Real alphaexp = 1.;
  Real betaexp = 2.;
  if (D < damage_limit_with_damage_rate) {
    Real alpha = std::max(0., 1. - pow((r * r / this->R2), alphaexp));
    w = pow(alpha, betaexp);
  }

  return w;
}

} // namespace akantu
