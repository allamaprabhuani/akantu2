/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_MAZARS_NON_LOCAL_TMPL_HH_
#define AKANTU_MATERIAL_MAZARS_NON_LOCAL_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
MaterialMazarsNonLocal<dim, Parent>::MaterialMazarsNonLocal(
    SolidMechanicsModel & model, const ID & id)
    : parent(model, id), Ehat(this->registerInternal("Ehat", 1)),
      non_local_variable(this->registerInternal("non_local_variable", 1)) {
  AKANTU_DEBUG_IN();

  this->is_non_local = true;

  this->registerParam("average_on_damage", this->damage_in_compute_stress,
                      false, _pat_parsable | _pat_modifiable,
                      "Is D the non local variable");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialMazarsNonLocal<dim, Parent>::registerNonLocalVariables() {
  ID local;
  if (this->damage_in_compute_stress) {
    local = this->damage.getName();
  } else {
    local = this->Ehat.getName();
  }

  this->getModel().getNonLocalManager().registerNonLocalVariable(
      local, non_local_variable.getName(), 1);
  this->getModel()
      .getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(non_local_variable.getName());
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialMazarsNonLocal<dim, Parent>::computeStress(ElementType el_type,
                                                        GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);

  for (auto && data : arguments) {
    parent::MaterialParent::computeStressOnQuad(data);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialMazarsNonLocal<dim, Parent>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {
  auto & non_loc_var = non_local_variable(el_type, ghost_type);

  if (this->damage_in_compute_stress) {
    auto && arguments = zip_replace(getArguments(el_type, ghost_type),
                                    "damage"_n = make_view(non_loc_var));

    for (auto && data : arguments) {
      parent::MaterialParent::computeDamageAndStressOnQuad(data);
    }
  } else {
    auto && arguments = zip_replace(getArguments(el_type, ghost_type),
                                    "Ehat"_n = make_view(non_loc_var));

    for (auto && data : arguments) {
      parent::MaterialParent::computeDamageAndStressOnQuad(data);
    }
  }
}

} // namespace akantu

#endif // AKANTU_MATERIAL_MAZARS_NON_LOCAL_TMPL_HH_
