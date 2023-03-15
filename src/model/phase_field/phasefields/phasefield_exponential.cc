/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "phasefield_exponential.hh"
#include "aka_common.hh"
#include <tuple>

namespace akantu {

/* -------------------------------------------------------------------------- */
PhaseFieldExponential::PhaseFieldExponential(PhaseFieldModel & model,
                                             const ID & id)
    : PhaseField(model, id) {}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::updateInternalParameters() {
  PhaseField::updateInternalParameters();

  for (const auto & type :
       element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    for (auto && tuple : zip(make_view(this->damage_energy(type, _not_ghost),
                                       spatial_dimension, spatial_dimension),
                             this->g_c(type, _not_ghost))) {
      Matrix<Real> d =
          Matrix<Real>::Identity(spatial_dimension, spatial_dimension) *
          std::get<1>(tuple) * this->l0;
      std::get<0>(tuple) = d;
    }
  }
}

/* -------------------------------------------------------------------------- */
void PhaseFieldExponential::computeDrivingForce(ElementType el_type,
                                                GhostType ghost_type) {
  for (auto && tuple : zip(this->phi(el_type, ghost_type),
                           this->phi.previous(el_type, ghost_type),
                           this->driving_force(el_type, ghost_type),
                           this->damage_energy_density(el_type, ghost_type),
                           make_view(this->strain(el_type, ghost_type),
                                     spatial_dimension, spatial_dimension),
                           this->g_c(el_type, ghost_type))) {
    computePhiOnQuad(std::get<4>(tuple), std::get<0>(tuple),
                     std::get<1>(tuple));
    computeDamageEnergyDensityOnQuad(std::get<0>(tuple), std::get<3>(tuple),
                                     std::get<5>(tuple));
    computeDrivingForceOnQuad(std::get<0>(tuple), std::get<2>(tuple));
  }
}

INSTANTIATE_PHASEFIELD(exponential, PhaseFieldExponential);

} // namespace akantu
