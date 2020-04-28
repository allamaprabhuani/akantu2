/**
 * @file   material_cohesive_linear_slipweakening.cc
 *
 * @author Mathias Lebihain <mathias.lebihain@epfl.ch>
 *
 * @date creation: Fri Mar 06 2020
 * @date last modification: Fri Mar 06 2020
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * Mohr-Coulomb insertion criterion and slip-weakening friction for
 * extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_friction_slipweakening.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialCohesiveLinearFrictionSlipWeakening<dim>::
    MaterialCohesiveLinearFrictionSlipWeakening(SolidMechanicsModel & model,
                                                const ID & id)
    : MaterialParent(model, id), mu_dynamic("mu_dynamic", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("mu_static", this->mu_insertion,
                      _pat_parsable | _pat_readable,
                      "Static value of the friction coefficient");

  this->registerParam("mu_dynamic", this->mu_dynamic,
                      _pat_parsable | _pat_readable,
                      "Dynamic value of the friction coefficient");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialCohesiveLinearFrictionSlipWeakening<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialParent::initMaterial();

  mu_dynamic.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/** \todo Frictional contributions might be computed from the contact_tractions,
 * and delta_max at the previous time step
 */
template <UInt dim>
void MaterialCohesiveLinearFrictionSlipWeakening<dim>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Vector<Real> normal_opening(dim);
  Vector<Real> tangential_opening(dim);

  /// loop on each quadrature point
  for (auto && data :
       zip(make_view(this->insertion_stress(el_type, ghost_type), dim),
           make_view(this->insertion_compression(el_type, ghost_type), dim),
           this->sigma_c_eff(el_type, ghost_type),
           this->delta_c_eff(el_type, ghost_type),
           this->mu_eff(el_type, ghost_type),
           this->mu_dynamic(el_type, ghost_type), make_view(normal, dim),
           this->delta_max(el_type, ghost_type),
           this->damage(el_type, ghost_type),
           make_view(this->tractions(el_type, ghost_type), dim),
           make_view(this->contact_tractions(el_type, ghost_type), dim),
           make_view(this->friction_force(el_type, ghost_type), dim),
           make_view(this->opening(el_type, ghost_type), dim),
           make_view(this->contact_opening(el_type, ghost_type), dim),
           this->residual_sliding(el_type, ghost_type),
           this->residual_sliding.previous(el_type, ghost_type))) {

    auto && insertion_traction = std::get<0>(data);
    auto && insertion_compression = std::get<1>(data);
    auto && sigma_c = std::get<2>(data);
    auto && delta_c = std::get<3>(data);
    auto && mu_static = std::get<4>(data);
    auto && mu_dynamic = std::get<5>(data);
    auto && normal = std::get<6>(data);
    auto && delta_max = std::get<7>(data);
    auto && damage = std::get<8>(data);
    auto && traction = std::get<9>(data);
    auto && contact_traction = std::get<10>(data);
    auto && friction_force = std::get<11>(data);
    auto && opening = std::get<12>(data);
    auto && contact_opening = std::get<13>(data);
    auto && res_sliding = std::get<14>(data);
    auto && res_sliding_prev = std::get<15>(data);

    Real normal_opening_norm, tangential_opening_norm;
    bool penetration;

    this->computeTractionOnQuad(
        traction, opening, normal, delta_max, delta_c, insertion_traction,
        insertion_compression, sigma_c, normal_opening, tangential_opening,
        normal_opening_norm, tangential_opening_norm, damage, penetration,
        contact_traction, contact_opening);

    if (penetration) {
      /// current friction coefficient
      Real mu = mu_static + (mu_dynamic - mu_static) * damage;
      /// the definition of tau_max refers to
      /// the current contact_traction computed in computeTractionOnQuad
      Real tau_max = mu * contact_traction.norm();

      // elastic sliding
      Real delta_sliding_norm =
          std::abs(tangential_opening_norm - res_sliding_prev);

      /// tau is the norm of the friction force, acting tangentially to the
      /// surface
      Real tau = std::min(this->friction_penalty * delta_sliding_norm, tau_max);

      if ((tangential_opening_norm - res_sliding_prev) < 0.0)
        tau = -tau;

      /// from tau get the x and y components of friction, to be added in the
      /// force vector
      auto tangent_unit_vector = tangential_opening / tangential_opening_norm;
      friction_force = tau * tangent_unit_vector;

      /// update residual_sliding
      res_sliding =
          tangential_opening_norm - (std::abs(tau) / this->friction_penalty);
    } else {
      friction_force.clear();
    }
    traction += friction_force;
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_linear_friction_slipweakening,
                     MaterialCohesiveLinearFrictionSlipWeakening);

} // namespace akantu
