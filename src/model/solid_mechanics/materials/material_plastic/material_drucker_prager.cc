/**
 * @file   material_drucker_prager.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Tue Apr 06 2021
 *
 * @brief  Implementation of the akantu::MaterialDruckerPrager class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_drucker_prager.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template <UInt spatial_dimension>
MaterialDruckerPrager<spatial_dimension>::MaterialDruckerPrager(
    SolidMechanicsModel & model, const ID & id)
    : MaterialPlastic<spatial_dimension>(model, id) {

  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialDruckerPrager<spatial_dimension>::MaterialDruckerPrager(
    SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
    FEEngine & fe_engine, const ID & id)
    : MaterialPlastic<spatial_dimension>(model, dim, mesh, fe_engine, id) {

  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDruckerPrager<spatial_dimension>::initialize() {
  this->registerParam("phi", phi, Real(0.), _pat_parsable | _pat_modifiable,
                      "Internal friction angle in degrees");
  this->registerParam("fc", fc, Real(1.), _pat_parsable | _pat_modifiable,
                      "Compressive strength");
  this->registerParam("radial_return", radial_return_mapping, bool(true),
                      _pat_parsable | _pat_modifiable, "Radial return mapping");

  this->updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDruckerPrager<spatial_dimension>::updateInternalParameters() {
  MaterialElastic<spatial_dimension>::updateInternalParameters();

  // compute alpha and k parameters for Drucker-Prager
  Real phi_radian = this->phi * M_PI / 180.;
  this->alpha = (6. * sin(phi_radian)) / (3. - sin(phi_radian));
  Real cohesion = this->fc * (1. - sin(phi_radian)) / (2. * cos(phi_radian));
  this->k = (6. * cohesion * cos(phi_radian)) / (3. - sin(phi_radian));
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDruckerPrager<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {

  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);
  // infinitesimal and finite deformation
  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_sigma_th_it =
      this->sigma_th.previous(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto previous_stress_it = this->stress.previous(el_type, ghost_type)
                                .begin(spatial_dimension, spatial_dimension);

  auto inelastic_strain_it = this->inelastic_strain(el_type, ghost_type)
                                 .begin(spatial_dimension, spatial_dimension);

  auto previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);

  //
  // Finite Deformations
  //
  if (this->finite_deformation) {
    auto previous_piola_kirchhoff_2_it =
        this->piola_kirchhoff_2.previous(el_type, ghost_type)
            .begin(spatial_dimension, spatial_dimension);

    auto green_strain_it = this->green_strain(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_piola_kirchhoff_2_it;

    auto & green_strain = *green_strain_it;
    this->template gradUToE<spatial_dimension>(grad_u, green_strain);
    Matrix<Real> previous_green_strain(spatial_dimension, spatial_dimension);
    this->template gradUToE<spatial_dimension>(previous_grad_u,
                                               previous_green_strain);
    Matrix<Real> F_tensor(spatial_dimension, spatial_dimension);
    this->template gradUToF<spatial_dimension>(grad_u, F_tensor);

    computeStressOnQuad(green_strain, previous_green_strain, sigma,
                        previous_sigma, inelastic_strain_tensor,
                        previous_inelastic_strain_tensor, *sigma_th_it,
                        *previous_sigma_th_it, F_tensor);

    ++sigma_th_it;
    ++inelastic_strain_it;
    ++previous_sigma_th_it;
    //++previous_stress_it;
    ++previous_gradu_it;
    ++green_strain_it;
    ++previous_inelastic_strain_it;
    ++previous_piola_kirchhoff_2_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  }
  // Infinitesimal deformations
  else {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_stress_it;

    computeStressOnQuad(
        grad_u, previous_grad_u, sigma, previous_sigma, inelastic_strain_tensor,
        previous_inelastic_strain_tensor, *sigma_th_it, *previous_sigma_th_it);
    ++sigma_th_it;
    ++inelastic_strain_it;
    ++previous_sigma_th_it;
    ++previous_stress_it;
    ++previous_gradu_it;
    ++previous_inelastic_strain_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDruckerPrager<spatial_dimension>::computeTangentModuli(
    ElementType /*el_type*/, Array<Real> & /*tangent_matrix*/,
    GhostType /*ghost_type*/) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(plastic_drucker_prager, MaterialDruckerPrager);

} // namespace akantu
