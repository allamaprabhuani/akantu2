/**
 * @file   material_phasefield.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Fri Apr 02 2021
 *
 * @brief  Specialization of the material class for the phasefield material
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
#include "material_phasefield.hh"
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialPhaseField<spatial_dimension>::MaterialPhaseField(
    SolidMechanicsModel & model, const ID & id)
    : Parent(model, id), effective_damage("effective_damage", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam("eta", eta, Real(0.), _pat_parsable, "eta");
  this->damage.initialize(0);
  this->effective_damage.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  computeEffectiveDamage(el_type, ghost_type);
  auto dam = this->effective_damage(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);
  sigma *= (1. - *dam) * (1. - *dam);

  ++dam;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeTangentModuli(el_type, tangent_matrix, ghost_type);

  computeEffectiveDamage(el_type, ghost_type);
  auto dam = this->effective_damage(el_type, ghost_type).begin();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  tangent *= (1. - *dam) * (1. - *dam) + eta;
  ++dam;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeEffectiveDamage(
    ElementType el_type, GhostType ghost_type) {

  auto && grad_u_view =
      make_view(this->gradu(el_type, ghost_type), this->spatial_dimension,
                this->spatial_dimension);

  for (auto && data : zip(grad_u_view, this->damage(el_type, ghost_type),
                          this->effective_damage(el_type, ghost_type))) {
    auto & grad_u = std::get<0>(data);
    auto & dam = std::get<1>(data);
    auto & eff_dam = std::get<2>(data);

    computeEffectiveDamageOnQuad(grad_u, dam, eff_dam);
  }
}

INSTANTIATE_MATERIAL(phasefield, MaterialPhaseField);

} // namespace akantu
