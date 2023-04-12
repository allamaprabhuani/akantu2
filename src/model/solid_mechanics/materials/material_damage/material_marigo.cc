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
#include "material_marigo.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialMarigo<dim>::MaterialMarigo(SolidMechanicsModel & model, const ID & id)
    : MaterialDamage<dim>(model, id) {
  this->Yd =
      this->template registerInternal<Real, RandomInternalField>("Yd", 1);

  this->registerParam("Sd", Sd, Real(5000.), _pat_parsable | _pat_modifiable);
  this->registerParam("epsilon_c", epsilon_c, Real(0.), _pat_parsable,
                      "Critical strain");
  this->registerParam("Yc limit", yc_limit, false, _pat_internal,
                      "As the material a critical Y");
  this->registerParam("damage_in_y", damage_in_y, false, _pat_parsable,
                      "Use threshold (1-D)Y");
  this->registerParam("Yd", *Yd, _pat_parsable, "Damaging energy threshold");
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialMarigo<dim>::initMaterial() {
  MaterialDamage<dim>::initMaterial();
  updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialMarigo<dim>::updateInternalParameters() {
  MaterialDamage<dim>::updateInternalParameters();
  Yc = epsilon_c * this->E * epsilon_c / 2.;
  yc_limit = (std::abs(epsilon_c) > std::numeric_limits<Real>::epsilon());
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMarigo<dim>::computeStress(ElementType el_type,
                                        GhostType ghost_type) {
  for (auto && data : getArguments(el_type, ghost_type)) {
    computeStressOnQuad(data);
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialMarigo<1>;
template class MaterialMarigo<2>;
template class MaterialMarigo<3>;

static bool material_is_alocated_marigo =
    instantiateMaterial<MaterialMarigo>("marigo");

} // namespace akantu
