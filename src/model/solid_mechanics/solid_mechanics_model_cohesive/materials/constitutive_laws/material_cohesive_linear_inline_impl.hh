/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include "material_cohesive_linear.hh"
#include "aka_static_if.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class D1, class D2, class D3, class D4>
Real MaterialCohesiveLinear<dim>::computeEffectiveNorm(
    const Eigen::MatrixBase<D1> & stress, const Eigen::MatrixBase<D2> & normal,
    const Eigen::MatrixBase<D3> & tangent,
    const Eigen::MatrixBase<D4> & normal_traction_) const {
  Eigen::MatrixBase<D4> & normal_traction =
      const_cast<Eigen::MatrixBase<D4> &>(normal_traction_);
  normal_traction = stress * normal;

  Real normal_contrib = normal_traction.dot(normal);

  /// in 3D tangential components must be summed
  Real tangent_contrib = 0;

  if constexpr (dim == 2) {
    Real tangent_contrib_tmp = normal_traction.dot(tangent);
    tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
  } else if constexpr (dim == 3) {
    for (auto && tangent_v : tangent) {
      Real tangent_contrib_tmp = normal_traction.dot(tangent_v);
      tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
    }
  }

  tangent_contrib = std::sqrt(tangent_contrib);
  normal_contrib = std::max(Real(0.), normal_contrib);

  return std::sqrt(normal_contrib * normal_contrib +
                   tangent_contrib * tangent_contrib * beta2_inv);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialCohesiveLinear<dim>::computeTractionOnQuad(Args && args) {
  auto && delta_c = args["delta_c"_n];
  auto && sigma_c = args["sigma_c"_n];
  auto && delta_max = args["delta_max"_n];

  auto && normal = args["normal"_n];
  auto && opening = args["opening"_n];
  auto && traction = args["traction"_n];

  auto && damage = args["damage"_n];
  auto && insertion_stress = args["insertion_stress"_n];

  /// compute normal and tangential opening vectors
  auto norm = opening.dot(normal);
  this->normal_opening_norm = norm;
  this->normal_opening = normal * this->normal_opening_norm;

  this->tangential_opening = opening - normal_opening;
  this->tangential_opening_norm = this->tangential_opening.norm();

  /**
   * compute effective opening displacement
   * @f$ \delta = \sqrt{
   * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
   */
  auto delta = this->tangential_opening_norm * this->tangential_opening_norm *
               this->beta2_kappa2;

  penetration = this->normal_opening_norm / delta_c < -Math::getTolerance();

  // penetration = normal_opening_norm < 0.;
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  auto && contact_traction = args.get("contact_traction"_n);
  auto && contact_opening = args.get("contact_opening"_n);

  if (penetration) {
    /// use penalty coefficient in case of penetration
    contact_traction = this->normal_opening * this->penalty;
    contact_opening = this->normal_opening;

    /// don't consider penetration contribution for delta
    opening = tangential_opening;
    this->normal_opening.zero();
  } else {
    delta += this->normal_opening_norm * this->normal_opening_norm;
    contact_traction.zero();
    contact_opening.zero();
  }

  delta = std::sqrt(delta);

  /// update maximum displacement and damage
  delta_max = std::max(delta_max, delta);
  damage = std::min(delta_max / delta_c, Real(1.));

  /**
   * Compute traction @f$ \mathbf{T} = \left(
   * \frac{\beta^2}{\kappa} \Delta_t \mathbf{t} + \Delta_n
   * \mathbf{n} \right) \frac{\sigma_c}{\delta} \left( 1-
   * \frac{\delta}{\delta_c} \right)@f$
   */

  if (Math::are_float_equal(damage, 1.)) {
    traction.zero();
  } else if (Math::are_float_equal(damage, 0.)) {
    if (penetration) {
      traction.zero();
    } else {
      traction = insertion_stress;
    }
  } else {
    AKANTU_DEBUG_ASSERT(delta_max != 0.,
                        "Division by zero, tolerance might be too low");

    traction =
        (this->tangential_opening * this->beta2_kappa + this->normal_opening) *
        sigma_c / delta_max * (1. - damage);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Derived, class Args>
inline void MaterialCohesiveLinear<dim>::computeTangentTractionOnQuad(
    Eigen::MatrixBase<Derived> & tangent, Args && args) {
  /**
   * During the update of the residual the interpenetrations are
   * stored in the array "contact_opening", therefore, in the case
   * of penetration, in the array "opening" there are only the
   * tangential components.
   */

  auto && delta_c = args["delta_c"_n];
  auto && delta_max = args["delta_max"_n];

  auto && normal = args["normal"_n];
  auto && opening = args["opening"_n];
  auto && contact_opening = args["contact_opening"_n];

  auto && damage = args["damage"_n];

  opening += contact_opening;

  /// compute normal and tangential opening vectors
  this->normal_opening_norm = opening.dot(normal);
  this->normal_opening = normal * normal_opening_norm;

  this->tangential_opening = opening - normal_opening;
  this->tangential_opening_norm = tangential_opening.norm();

  auto delta = this->tangential_opening_norm * this->tangential_opening_norm *
               this->beta2_kappa2;

  penetration = normal_opening_norm < 0.0;
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  Real derivative = 0; // derivative = d(t/delta)/ddelta
  Real t = 0;

  auto && n_outer_n = normal * normal.transpose();

  if (penetration) {
    /// stiffness in compression given by the penalty parameter
    tangent = penalty * (tangent + n_outer_n);

    opening = tangential_opening;

    normal_opening_norm = opening.dot(normal);
    normal_opening = normal * normal_opening_norm;
  } else {
    delta += normal_opening_norm * normal_opening_norm;
  }

  delta = std::sqrt(delta);

  /**
   * Delta has to be different from 0 to have finite values of
   * tangential stiffness.  At the element insertion, delta =
   * 0. Therefore, a fictictious value is defined, for the
   * evaluation of the first value of K.
   */
  if (delta < Math::getTolerance()) {
    delta = delta_c / 1000.;
  }

  if (delta >= delta_max) {
    if (delta <= delta_c) {
      derivative = -sigma_c / (delta * delta);
      t = sigma_c * (1 - delta / delta_c);
    } else {
      derivative = 0.;
      t = 0.;
    }
  } else if (delta < delta_max) {
    Real tmax = sigma_c * (1 - delta_max / delta_c);
    t = tmax / delta_max * delta;
  }

  /// computation of the derivative of the constitutive law (dT/ddelta)
  auto && I = Eigen::Matrix<Real, dim, dim>::Identity() * this->beta2_kappa;

  auto && nn = (n_outer_n * (1. - this->beta2_kappa) + I) * t / delta;
  auto && mm = opening * this->beta2_kappa2;

  auto && t_tilde = normal_opening * (1. - this->beta2_kappa2) + mm;
  auto && t_hat = normal_opening + this->beta2_kappa * tangential_opening;

  auto && prov = t_hat * t_tilde.transpose() * derivative / delta + nn;

  tangent += prov.transpose();
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif // AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_HH_
