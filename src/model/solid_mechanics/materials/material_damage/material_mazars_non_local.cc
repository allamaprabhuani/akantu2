/**
 * @file   material_mazars_non_local.cc
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Specialization of the material class for the non-local mazars
 * material
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
#include "material_mazars_non_local.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialMazarsNonLocal<dim>::MaterialMazarsNonLocal(SolidMechanicsModel & model,
                                                    const ID & id)
    : parent(model, id), Ehat("epsilon_equ", *this),
      non_local_variable("mazars_non_local", *this) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;
  this->Ehat.initialize(1);
  this->non_local_variable.initialize(1);

  this->registerParam("average_on_damage", this->damage_in_compute_stress,
                      false, _pat_parsable | _pat_modifiable,
                      "Is D the non local variable");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMazarsNonLocal<dim>::registerNonLocalVariables() {
  ID local;
  if (this->damage_in_compute_stress) {
    local = this->damage.getName();
  } else {
    local = this->Ehat.getName();
  }

  this->model.getNonLocalManager().registerNonLocalVariable(
      local, non_local_variable.getName(), 1);
  this->model.getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(non_local_variable.getName());
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMazarsNonLocal<dim>::computeStress(ElementType el_type,
                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto && arguments = getArguments(el_type, ghost_type);

  for (auto && data : arguments) {
    MaterialMazars<dim>::computeStressOnQuad(data);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMazarsNonLocal<dim>::computeNonLocalStress(ElementType el_type,
                                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & non_loc_var = non_local_variable(el_type, ghost_type);

  if (this->damage_in_compute_stress) {
    auto && arguments = zip_replace<"damage"_h>(
        getArguments(el_type, ghost_type), make_view(non_loc_var));

    for (auto && data : arguments) {
      MaterialMazars<dim>::computeDamageAndStressOnQuad(data);
    }
  } else {
    auto && arguments = zip_replace<"Ehat"_h>(getArguments(el_type, ghost_type),
                                              make_view(non_loc_var));

    for (auto && data : arguments) {
      MaterialMazars<dim>::computeDamageAndStressOnQuad(data);
    }
  }
  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(mazars_non_local, MaterialMazarsNonLocal);

} // namespace akantu
