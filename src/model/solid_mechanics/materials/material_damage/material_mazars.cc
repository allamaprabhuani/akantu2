/**
 * @file   material_mazars.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Jul 24 2020
 *
 * @brief  Specialization of the material class for the damage material
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
#include "material_mazars.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialMazars<dim>::MaterialMazars(SolidMechanicsModel & model, const ID & id)
    : MaterialDamage<dim>(model, id), K0("K0", *this),
      damage_in_compute_stress(true) {
  AKANTU_DEBUG_IN();

  this->registerParam("K0", K0, _pat_parsable, "K0");
  this->registerParam("At", At, Real(0.8), _pat_parsable, "At");
  this->registerParam("Ac", Ac, Real(1.4), _pat_parsable, "Ac");
  this->registerParam("Bc", Bc, Real(1900.), _pat_parsable, "Bc");
  this->registerParam("Bt", Bt, Real(12000.), _pat_parsable, "Bt");
  this->registerParam("beta", beta, Real(1.06), _pat_parsable, "beta");

  this->K0.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialMazars<dim>::computeStress(ElementType el_type,
                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & damage = this->damage(el_type, ghost_type);
  Real Ehat = 0;

  auto && arguments = named_zip(
      tuple::get<"damage"_h>() = make_view(damage),
      tuple::get<"sigma"_h>() = make_view<dim, dim>(this->stress(el_type, ghost_type)),
      tuple::get<"grad_u"_h>() = make_view<dim, dim>(this->grad_u(el_type, ghost_type)),
      tuple::get<"Ehat"_h>() = broadcast(Ehat, damage.size()));

  for(auto && data : arguments) {
    computeStressOnQuad(data);
  }
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(mazars, MaterialMazars);

} // namespace akantu
