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
#include "material_marigo_non_local.hh"
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialMarigoNonLocal<dim>::MaterialMarigoNonLocal(SolidMechanicsModel & model,
                                                    const ID & id)
    : parent(model, id), Y(this->template registerInternal<Real>("Y", 1)),
      Ynl(this->template registerInternal<Real>("Y non local", 1)) {
  this->is_non_local = true;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMarigoNonLocal<dim>::registerNonLocalVariables() {
  this->getModel().getNonLocalManager().registerNonLocalVariable(
      this->Y.getName(), Ynl.getName(), 1);
  this->getModel()
      .getNonLocalManager()
      .getNeighborhood(this->name)
      .registerNonLocalVariable(Ynl.getName());
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMarigoNonLocal<dim>::computeStress(ElementType el_type,
                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto && arguments = getArguments(el_type, ghost_type);

  for (auto && data : arguments) {
    MaterialMarigo<dim>::computeStressOnQuad(data);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMarigoNonLocal<dim>::computeNonLocalStress(ElementType el_type,
                                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto && arguments = getArgumentsNonLocal(el_type, ghost_type);

  for (auto && data : arguments) {
    this->computeDamageAndStressOnQuad(data);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class MaterialMarigoNonLocal<1>;
template class MaterialMarigoNonLocal<2>;
template class MaterialMarigoNonLocal<3>;

const bool material_is_alocated_marigo_non_local [[maybe_unused]] =
    instantiateMaterial<MaterialMarigo>("marigo_non_local");

} // namespace akantu
