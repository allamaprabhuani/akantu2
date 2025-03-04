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
#include <algorithm>
#include <numeric>

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_uncoupled.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveLinearUncoupled<dim>::MaterialCohesiveLinearUncoupled(
    SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<dim>(model, id),
      delta_n_max(this->template registerInternal<Real, CohesiveInternalField>(
          "delta_n_max", 1)),
      delta_t_max(this->template registerInternal<Real, CohesiveInternalField>(
          "delta_t_max", 1)),
      damage_n(this->template registerInternal<Real, CohesiveInternalField>(
          "damage_n", 1)),
      damage_t(this->template registerInternal<Real, CohesiveInternalField>(
          "damage_t", 1)) {

  AKANTU_DEBUG_IN();

  this->registerParam(
      "roughness", R, Real(1.), _pat_parsable | _pat_readable,
      "Roughness to define coupling between mode II and mode I");

  delta_n_max.initializeHistory();
  delta_t_max.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinearUncoupled<dim>::computeTraction(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  delta_n_max.resize();
  delta_t_max.resize();
  damage_n.resize();
  damage_t.resize();

  /// define iterators
  auto traction_it = this->tractions(el_type, ghost_type).begin(dim);

  auto traction_end = this->tractions(el_type, ghost_type).end(dim);

  auto opening_it = this->opening(el_type, ghost_type).begin(dim);
  auto contact_traction_it =
      this->contact_tractions(el_type, ghost_type).begin(dim);
  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(dim);

  auto normal_it = this->normals(el_type, ghost_type).begin(dim);
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_n_max_it = delta_n_max(el_type, ghost_type).begin();
  auto delta_t_max_it = delta_t_max(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_n_it = damage_n(el_type, ghost_type).begin();
  auto damage_t_it = damage_t(el_type, ghost_type).begin();

  auto insertion_stress_it =
      this->insertion_stress(el_type, ghost_type).begin(dim);

  Vector<Real, dim> normal_opening;
  Vector<Real, dim> tangential_opening;

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++contact_traction_it, ++contact_opening_it,
       ++normal_it, ++sigma_c_it, ++delta_n_max_it, ++delta_t_max_it,
       ++delta_c_it, ++damage_n_it, ++damage_t_it, ++insertion_stress_it) {

    Real normal_opening_norm;
    Real tangential_opening_norm;
    bool penetration;

    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / R / R;

    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening = *normal_it * normal_opening_norm;

    tangential_opening = *opening_it - normal_opening;
    tangential_opening_norm = tangential_opening.norm();

    /// compute effective opening displacement
    Real delta_n =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    if (not this->contact_after_breaking and
        Math::are_float_equal(*damage_n_it, 1.)) {
      penetration = false;
    }

    if (penetration) {
      /// use penalty coefficient in case of penetration
      *contact_traction_it = normal_opening * this->penalty;
      *contact_opening_it = normal_opening;

      /// don't consider penetration contribution for delta
      //*opening_it = tangential_opening;
      normal_opening.zero();

    } else {
      delta_n += normal_opening_norm * normal_opening_norm;
      delta_t += normal_opening_norm * normal_opening_norm * delta_c2_R2;
      contact_traction_it->zero();
      contact_opening_it->zero();
    }

    delta_n = std::sqrt(delta_n);
    delta_t = std::sqrt(delta_t);

    /// update maximum displacement and damage
    *delta_n_max_it = std::max(*delta_n_max_it, delta_n);
    *damage_n_it = std::min(*delta_n_max_it / *delta_c_it, Real(1.));

    *delta_t_max_it = std::max(*delta_t_max_it, delta_t);
    *damage_t_it = std::min(*delta_t_max_it / *delta_c_it, Real(1.));

    Vector<Real, dim> normal_traction;
    Vector<Real, dim> shear_traction;

    /// NORMAL TRACTIONS
    if (Math::are_float_equal(*damage_n_it, 1.)) {
      normal_traction.zero();
    } else if (Math::are_float_equal(*damage_n_it, 0.)) {
      if (penetration) {
        normal_traction.zero();
      } else {
        normal_traction = *insertion_stress_it;
      }
    } else {
      // the following formulation holds both in loading and in
      // unloading-reloading
      AKANTU_DEBUG_ASSERT(*delta_n_max_it != 0.,
                          "Division by zero, tolerance might be too low");

      normal_traction = normal_opening * *sigma_c_it / (*delta_n_max_it) *
                        (1. - *damage_n_it);
    }

    /// SHEAR TRACTIONS
    if (Math::are_float_equal(*damage_t_it, 1.) or
        Math::are_float_equal(*damage_t_it, 0.)) {
      shear_traction.zero();
    } else {
      AKANTU_DEBUG_ASSERT(*delta_t_max_it != 0.,
                          "Division by zero, tolerance might be too low");

      shear_traction = tangential_opening * this->beta2_kappa * *sigma_c_it /
                       (*delta_t_max_it) * (1. - *damage_t_it);
    }

    *traction_it = normal_traction + shear_traction;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinearUncoupled<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  auto tangent_it = tangent_matrix.begin(dim, dim);
  auto tangent_end = tangent_matrix.end(dim, dim);
  auto normal_it = this->normals(el_type, ghost_type).begin(dim);
  auto opening_it = this->opening(el_type, ghost_type).begin(dim);

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  auto delta_n_max_it = delta_n_max.previous(el_type, ghost_type).begin();
  auto delta_t_max_it = delta_t_max.previous(el_type, ghost_type).begin();
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_n_it = damage_n(el_type, ghost_type).begin();

  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(dim);

  Vector<Real, dim> normal_opening;
  Vector<Real, dim> tangential_opening;

  for (; tangent_it != tangent_end; ++tangent_it, ++normal_it, ++opening_it,
                                    ++sigma_c_it, ++delta_c_it,
                                    ++delta_n_max_it, ++delta_t_max_it,
                                    ++damage_n_it, ++contact_opening_it) {

    Real normal_opening_norm;
    Real tangential_opening_norm;
    bool penetration;
    Real delta_c2_R2 = *delta_c_it * (*delta_c_it) / R / R;

    /**
     * During the update of the residual the interpenetrations are
     * stored in the array "contact_opening", therefore, in the case
     * of penetration, in the array "opening" there are only the
     * tangential components.
     */
    *opening_it += *contact_opening_it;

    /// compute normal and tangential opening vectors
    normal_opening_norm = opening_it->dot(*normal_it);
    normal_opening = *normal_it * normal_opening_norm;

    tangential_opening = *opening_it - normal_opening;
    tangential_opening_norm = tangential_opening.norm();

    Real delta_n =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;
    Real delta_t =
        tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

    penetration = normal_opening_norm < 0.0;
    if (not this->contact_after_breaking and
        Math::are_float_equal(*damage_n_it, 1.)) {
      penetration = false;
    }

    Real derivative = 0; // derivative = d(t/delta)/ddelta
    Real T = 0;

    /// TANGENT STIFFNESS FOR NORMAL TRACTIONS
    Matrix<Real, dim, dim> n_outer_n = *normal_it * (*normal_it).transpose();

    if (penetration) {
      /// stiffness in compression given by the penalty parameter
      *tangent_it = n_outer_n * this->penalty;

      //*opening_it = tangential_opening;
      normal_opening.zero();
    } else {
      delta_n += normal_opening_norm * normal_opening_norm;
      delta_n = std::sqrt(delta_n);

      delta_t += normal_opening_norm * normal_opening_norm * delta_c2_R2;

      /**
       * Delta has to be different from 0 to have finite values of
       * tangential stiffness.  At the element insertion, delta =
       * 0. Therefore, a fictictious value is defined, for the
       * evaluation of the first value of K.
       */
      if (delta_n < Math::getTolerance()) {
        delta_n = *delta_c_it / 1000.;
      }

      // loading
      if (delta_n >= *delta_n_max_it) {
        if (delta_n <= *delta_c_it) {
          derivative = -(*sigma_c_it) / (delta_n * delta_n);
          T = *sigma_c_it * (1 - delta_n / *delta_c_it);
        } else {
          derivative = 0.;
          T = 0.;
        }
        // unloading-reloading
      } else if (delta_n < *delta_n_max_it) {
        Real T_max = *sigma_c_it * (1 - *delta_n_max_it / *delta_c_it);
        derivative = 0.;
        T = T_max / *delta_n_max_it * delta_n;
      }

      /// computation of the derivative of the constitutive law (dT/ddelta)
      Matrix<Real, dim, dim> nn = n_outer_n * T / delta_n;

      Vector<Real, dim> Delta_tilde =
          normal_opening * (1. - this->beta2_kappa2) +
          *opening_it * this->beta2_kappa2;

      Matrix<Real, dim, dim> prov =
          normal_opening * Delta_tilde.transpose() * derivative / delta_n + nn;

      *tangent_it = prov.transpose();
    }

    derivative = 0.;
    T = 0.;

    /// TANGENT STIFFNESS FOR SHEAR TRACTIONS
    delta_t = std::sqrt(delta_t);

    /**
     * Delta has to be different from 0 to have finite values of
     * tangential stiffness.  At the element insertion, delta =
     * 0. Therefore, a fictictious value is defined, for the
     * evaluation of the first value of K.
     */
    if (delta_t < Math::getTolerance()) {
      delta_t = *delta_c_it / 1000.;
    }

    // loading
    if (delta_t >= *delta_t_max_it) {
      if (delta_t <= *delta_c_it) {
        derivative = -(*sigma_c_it) / (delta_t * delta_t);
        T = *sigma_c_it * (1 - delta_t / *delta_c_it);
      } else {
        derivative = 0.;
        T = 0.;
      }
      // unloading-reloading
    } else if (delta_t < *delta_t_max_it) {
      Real T_max = *sigma_c_it * (1 - *delta_t_max_it / *delta_c_it);
      derivative = 0.;
      T = T_max / *delta_t_max_it * delta_t;
    }

    /// computation of the derivative of the constitutive law (dT/ddelta)
    auto I = (Matrix<Real, dim, dim>::Identity() - n_outer_n) * (T / delta_t);

    auto && Delta_tilde = normal_opening * (delta_c2_R2 - this->beta2_kappa2) +
                          *opening_it * this->beta2_kappa2;

    auto && Delta_hat = tangential_opening * this->beta2_kappa;

    auto && prov =
        Delta_hat * Delta_tilde.transpose() * derivative / delta_t + I;

    *tangent_it += prov.transpose();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveLinearUncoupled<1>;
template class MaterialCohesiveLinearUncoupled<2>;
template class MaterialCohesiveLinearUncoupled<3>;
const bool material_is_allocated_cohesive_linear_uncoupled [[maybe_unused]] =
    instantiateMaterial<MaterialCohesiveLinearUncoupled>(
        "cohesive_linear_uncoupled");

} // namespace akantu
