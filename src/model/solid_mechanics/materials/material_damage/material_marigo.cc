/**
 * @file   material_marigo.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Specialization of the material class for the marigo material
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
#include "material_marigo.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialMarigo<dim>::MaterialMarigo(SolidMechanicsModel & model, const ID & id)
    : MaterialDamage<dim>(model, id), Yd("Yd", *this), damage_in_y(false),
      yc_limit(false) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sd", Sd, Real(5000.), _pat_parsable | _pat_modifiable);
  this->registerParam("epsilon_c", epsilon_c, Real(0.), _pat_parsable,
                      "Critical strain");
  this->registerParam("Yc limit", yc_limit, false, _pat_internal,
                      "As the material a critical Y");
  this->registerParam("damage_in_y", damage_in_y, false, _pat_parsable,
                      "Use threshold (1-D)Y");
  this->registerParam("Yd", Yd, _pat_parsable, "Damaging energy threshold");

  this->Yd.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialMarigo<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialDamage<dim>::initMaterial();

  updateInternalParameters();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialMarigo<dim>::updateInternalParameters() {
  MaterialDamage<dim>::updateInternalParameters();
  Yc = .5 * epsilon_c * this->E * epsilon_c;
  yc_limit = (std::abs(epsilon_c) > std::numeric_limits<Real>::epsilon());
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMarigo<dim>::computeStress(ElementType el_type,
                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Y = 0.;

  auto && arguments = getArguments(el_type, ghost_type);

  for (auto && data : arguments) {
    computeStressOnQuad(data);
  }

  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(marigo, MaterialMarigo);

} // namespace akantu
