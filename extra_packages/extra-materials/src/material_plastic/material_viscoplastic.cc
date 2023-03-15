/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_viscoplastic.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialViscoPlastic<dim>::MaterialViscoPlastic(SolidMechanicsModel & model,
                                                const ID & id)
    : MaterialPlastic<dim>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("rate", rate, 0., _pat_parsable | _pat_modifiable,
                      "Rate sensitivity component");
  this->registerParam("edot0", edot0, 0., _pat_parsable | _pat_modifiable,
                      "Reference strain rate");
  this->registerParam("ts", ts, 0., _pat_parsable | _pat_modifiable,
                      "Time Step");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoPlastic<dim>::computeStress(ElementType el_type,
                                              GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real * iso_hardening = this->iso_hardening(el_type, ghost_type).data();

  auto previous_grad_u_it =
      this->gradu.previous(el_type, ghost_type).begin(dim, dim);

  auto previous_sigma_it =
      this->stress.previous(el_type, ghost_type).begin(dim, dim);

  auto inelastic_strain_it =
      this->inelastic_strain(el_type, ghost_type).begin(dim, dim);

  auto previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type, ghost_type).begin(dim, dim);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeStressOnQuad(grad_u, *previous_grad_u_it, sigma, *previous_sigma_it,
                      *inelastic_strain_it, *previous_inelastic_strain_it,
                      *iso_hardening);

  ++inelastic_strain_it;
  ++iso_hardening;
  ++previous_grad_u_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialViscoPlastic<dim>::computeTangentModuli(
    __attribute__((unused)) ElementType el_type, Array<Real> & tangent_matrix,
    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto previous_sigma_it =
      this->stress.previous(el_type, ghost_type).begin(dim, dim);

  auto previous_strain_it =
      this->gradu.previous(el_type, ghost_type).begin(dim, dim);

  Real * iso_hardening = this->iso_hardening(el_type, ghost_type).data();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  Matrix<Real> & previous_grad_u = *previous_strain_it;
  Matrix<Real> & previous_sigma_tensor = *previous_sigma_it;

  computeTangentModuliOnQuad(tangent, grad_u, previous_grad_u, sigma,
                             previous_sigma_tensor, *iso_hardening);
  ++previous_sigma_it;
  ++previous_strain_it;
  ++iso_hardening;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(visco_plastic, MaterialViscoPlastic);

} // namespace akantu
