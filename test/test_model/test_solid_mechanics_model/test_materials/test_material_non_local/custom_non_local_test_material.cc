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
#include "custom_non_local_test_material.hh"
#include "aka_types.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
CustomNonLocalTestMaterial<dim>::CustomNonLocalTestMaterial(
    SolidMechanicsModel & model, const ID & id)
    : MyNonLocalParent(model, id),
      local_damage(this->registerInternal("local_damage", 1)),
      damage(this->registerInternal("damage", 1)) {}

/* -------------------------------------------------------------------------- */
template <Int dim>
void CustomNonLocalTestMaterial<dim>::registerNonLocalVariables() {
  /// register the non-local variable in the manager
  this->getModel().getNonLocalManager().registerNonLocalVariable(
      this->local_damage.getName(), this->damage.getName(), 1);

  this->getModel()
      .getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(damage.getName());
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void CustomNonLocalTestMaterial<dim>::computeStress(ElementType el_type,
                                                    GhostType ghost_type) {
  MyNonLocalParent::computeStress(el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void CustomNonLocalTestMaterial<dim>::computeNonLocalStress(
    ElementType el_type, GhostType ghost_type) {

  // compute the damage and update the stresses
  for (auto && [stress, dam] :
       zip(make_view<dim, dim>(this->stress(el_type, ghost_type)),
           this->damage(el_type, ghost_type))) {
    stress = stress * (1. - dam);
  }
}

/* -------------------------------------------------------------------------- */
// Instantiate the material for the 3 dimensions
template class CustomNonLocalTestMaterial<1>;
template class CustomNonLocalTestMaterial<2>;
template class CustomNonLocalTestMaterial<3>;

const bool material_is_allocated_custom_non_local_test_material
    [[maybe_unused]] = instantiateMaterial<CustomNonLocalTestMaterial>(
        "custom_non_local_test_material");
/* -------------------------------------------------------------------------- */
} // namespace akantu
