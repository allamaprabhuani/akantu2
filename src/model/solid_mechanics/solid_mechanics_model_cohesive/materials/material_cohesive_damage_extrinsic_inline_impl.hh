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
#include "material_cohesive_damage_extrinsic.hh"

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialCohesiveDamageExtrinsic<dim>::computeTractionOnQuad(Args && args) {
  auto && d = args["czm_damage"_n];
  auto && opening = args["opening"_n];
  auto && traction = args["traction"_n];

  Real A = this->augmented_stiffness(d);
  traction = k*A*opening;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif // AKANTU_MATERIAL_COHESIVE_DAMAGE_EXTRINSIC_INLINE_IMPL_HH_
