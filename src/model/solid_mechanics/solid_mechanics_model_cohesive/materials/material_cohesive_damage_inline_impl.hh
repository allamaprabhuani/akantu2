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
//#include "material_cohesive_damage.hh"
#include "aka_static_if.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_DAMAGE_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_COHESIVE_DAMAGE_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialCohesiveDamage<dim>::computeTractionOnQuad(Args && args) {
    auto && lambda = args["lambda"_n];
    auto && err_opening = args["err_opening"_n];
    auto && opening = args["opening"_n];
    auto && traction = args["traction"_n];

    traction = lambda - (opening*k);
    /// TODO : COMPUTE augmented_compliance
    Real augmented_compliance(0.);
    err_opening = opening-lambda*augmented_compliance;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif // AKANTU_MATERIAL_COHESIVE_DAMAGE_INLINE_IMPL_HH_
