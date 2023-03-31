/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear_friction.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveLinearFriction<dim>::MaterialCohesiveLinearFriction(
    SolidMechanicsModel & model, const ID & id)
    : MaterialParent(model, id), residual_sliding("residual_sliding", *this),
      friction_force("friction_force", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("mu", mu_max, Real(0.), _pat_parsable | _pat_readable,
                      "Maximum value of the friction coefficient");

  this->registerParam("penalty_for_friction", friction_penalty, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Penalty parameter for the friction behavior");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveLinearFriction<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialParent::initMaterial();

  friction_force.initialize(dim);
  residual_sliding.initialize(1);
  residual_sliding.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
void MaterialCohesiveLinearFriction<dim>::computeTractionOnQuad(Args && args) {
  MaterialParent::computeTractionOnQuad(args);
  using namespace tuple;

  auto && res_sliding_prev = args["previous_residual_sliding"_n];
  auto && res_sliding = args["residual_sliding"_n];

  if (not this->penetration) {
    args["friction_force"_n].zero();
    return;
  }

  /// the friction coefficient mu increases with the damage. It
  /// equals the maximum value when damage = 1.
  //      Real damage = std::min(*delta_max_prev_it / *delta_c_it,
  //      Real(1.));
  Real mu = mu_max; // * damage;

  /// the definition of tau_max refers to the opening
  /// (penetration) of the previous incremental step

  Real tau_max = mu * this->penalty * (std::abs(this->normal_opening_norm));
  Real delta_sliding_norm =
      std::abs(this->tangential_opening_norm - res_sliding_prev);

  /// tau is the norm of the friction force, acting tangentially to the
  /// surface
  Real tau = std::min(friction_penalty * delta_sliding_norm, tau_max);

  if ((this->tangential_opening_norm - res_sliding_prev) < 0.0) {
    tau = -tau;
  }

  /// from tau get the x and y components of friction, to be added in the
  /// force vector
  args["friction_force"_n] =
      tau * this->tangential_opening / this->tangential_opening_norm;

  /// update residual_sliding
  if (friction_penalty == 0.) {
    res_sliding = this->tangential_opening_norm;
  } else {
    res_sliding =
        this->tangential_opening_norm - (std::abs(tau) / friction_penalty);
  }

  args["traction"_n] += args["friction_force"_n];
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinearFriction<dim>::computeTraction(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  using namespace tuple;

  if (not this->model->isDefaultSolverExplicit()) {
    this->delta_max(el_type, ghost_type)
        .copy(this->delta_max.previous(el_type, ghost_type));
  }

  for (auto && args : getArguments(el_type, ghost_type)) {
    computeTractionOnQuad(args);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinearFriction<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  for (auto && [args, tangent] : zip(getArguments(el_type, ghost_type),
                                     make_view<dim, dim>(tangent_matrix))) {
    this->computeTangentTractionOnQuad(tangent, args);

    if (not this->penetration) {
      continue;
    }
    //      Real damage = std::min(*delta_max_it / *delta_c_it, Real(1.));
    Real mu = mu_max; // * damage;

    Real previous_normal_opening_norm =
        std::min(args["previous_opening"_n].dot(args["normal"_n]), Real(0.));
    //      Vector<Real> normal_opening_prev = (*normal_it);
    //      normal_opening_prev *= normal_opening_prev_norm;

    Real tau_max =
        mu * this->penalty * (std::abs(previous_normal_opening_norm));
    Real delta_sliding_norm = std::abs(this->tangential_opening_norm -
                                       args["previous_residual_sliding"_n]);

    // tau is the norm of the friction force, acting tangentially to the
    // surface
    Real tau = std::min(this->friction_penalty * delta_sliding_norm, tau_max);

    if (tau < tau_max && tau_max > Math::getTolerance()) {
      auto && I = Matrix<Real, dim, dim>::Identity();
      auto && n = args["normal"_n];
      tangent += (I - n * n.transpose()) * friction_penalty;
    }
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template class MaterialCohesiveLinearFriction<1>;
template class MaterialCohesiveLinearFriction<2>;
template class MaterialCohesiveLinearFriction<3>;
static bool material_is_allocated_cohesive_linear_friction =
    instantiateMaterial<MaterialCohesiveLinearFriction>(
        "cohesive_linear_friction");

} // namespace akantu
