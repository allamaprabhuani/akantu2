/**
 * @file   material_cohesive_linear_inline_impl.hh
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Apr 22 2015
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Inline functions of the MaterialCohesiveLinear
 *
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
#include "material_cohesive_linear_sequential.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_COHESIVE_LINEAR_SEQUENTIAL_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/* ------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialCohesiveLinearSequential<dim>::computeTractionOnQuad(
    Vector<Real> & traction, Vector<Real> & opening,
    const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
    const Vector<Real> & insertion_stress, const Real & sigma_c,
    Vector<Real> & normal_opening, Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, Vector<Real> & contact_traction,
    Vector<Real> & contact_opening) {

  /// compute normal and tangential opening vectors
  normal_opening_norm = opening.dot(normal);
  normal_opening = normal;
  normal_opening *= normal_opening_norm;

  tangential_opening = opening;
  tangential_opening -= normal_opening;
  tangential_opening_norm = tangential_opening.norm();

  /**
   * compute effective opening displacement
   * @f$ \delta = \sqrt{
   * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
   */
  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  damage = std::min(delta_max / delta_c, Real(1.));
  penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  // penetration = normal_opening_norm < 0.;
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  if (penetration) {
    /// use penalty coefficient in case of penetration
    contact_traction = normal_opening;
    contact_traction *= this->penalty;
    contact_opening = normal_opening;

    /// don't consider penetration contribution for delta
    opening = tangential_opening;
    normal_opening.zero();
  } else {
    delta += normal_opening_norm * normal_opening_norm;
    contact_traction.zero();
    contact_opening.zero();
  }

  delta = std::sqrt(delta);

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
    traction = tangential_opening;
    traction *= this->beta2_kappa;
    traction += normal_opening;

    AKANTU_DEBUG_ASSERT(delta_max != 0.,
                        "Division by zero, tolerance might be too low");

    traction *= sigma_c / delta_max * (1. - damage);
  }
}

/* ------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialCohesiveLinearSequential<dim>::computeTangentTractionOnQuad(
    Matrix<Real> & tangent, Real & delta_max, const Real & delta_c,
    const Real & sigma_c, const Vector<Real> & normal,
    const Real & normal_opening_norm, const Real & tangential_opening_norm,
    const Real & damage) {

  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  auto penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  Real t = 0;

  Matrix<Real> n_outer_n(dim, dim);
  n_outer_n.outerProduct(normal, normal);

  if (penetration) {
    /// stiffness in compression given by the penalty parameter
    tangent += n_outer_n;
    tangent *= this->penalty;
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
  if (delta_max < Math::getTolerance()) {
    delta_max = delta_c / 100.;
  }
  if (delta < Math::getTolerance()) {
    delta = delta_c / 1000.;
  }

  if (delta_max <= delta_c) {
    Real tmax = sigma_c * (1 - delta_max / delta_c);
    t = tmax / delta_max;
  } else {
    t = 0.;
  }

  /// computation of the derivative of the constitutive law (dT/ddelta)
  Matrix<Real> I(dim, dim);
  I.eye(this->beta2_kappa);

  Matrix<Real> nn(n_outer_n);
  nn *= (1. - this->beta2_kappa);
  nn += I;
  nn *= t;

  Matrix<Real> prov = nn;
  Matrix<Real> prov_t = prov.transpose();
  if (not penetration)
    tangent += prov_t;
}

/* ------------------------------------------------------------------- */
template <UInt dim>
inline bool MaterialCohesiveLinearSequential<dim>::updateDeltaMaxOnQuad(
    const Real & normal_opening_norm, const Real & tangential_opening_norm,
    Real & damage, Real & delta_max, const Real & delta_c) {

  bool delta_max_increased{false};

  /// if damage is maximum, no point to increase delta_max
  if (Math::are_float_equal(damage, 1.))
    return delta_max_increased;

  /**
   * compute effective opening displacement
   * @f$ \delta = \sqrt{
   * \frac{\beta^2}{\kappa^2} \Delta_t^2 + \Delta_n^2 } @f$
   */
  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  auto penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  // penetration = normal_opening_norm < 0.;
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  if (not penetration) {
    delta += normal_opening_norm * normal_opening_norm;
  }

  delta = std::sqrt(delta);

  AKANTU_DEBUG_ASSERT(delta < (1 + 1e-3) * delta_max,
                      "Delta " << delta << " is way larger than delta_max "
                               << delta_max << " damage "
                               << delta_max / delta_c);

  /// update maximum displacement and damage
  if (std::abs(delta_max - delta) < 1e-3 * delta_max) {
    delta_max += delta_c / 10;
    delta_max_increased = true;
  }
  damage = std::min(delta_max / delta_c, Real(1.));

  return delta_max_increased;
}
/* ------------------------------------------------------------------- */
template <UInt dim>
inline Real MaterialCohesiveLinearSequential<dim>::computeDeltaMaxExcessOnQuad(
    const Real & normal_opening_norm, const Real & tangential_opening_norm,
    const Real & damage, const Real & delta_max, const Real & delta_c) {

  Real delta_max_excess(0);

  Real delta =
      tangential_opening_norm * tangential_opening_norm * this->beta2_kappa2;

  auto penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  if (not this->contact_after_breaking and Math::are_float_equal(damage, 1.)) {
    penetration = false;
  }

  if (not penetration) {
    delta += normal_opening_norm * normal_opening_norm;
  }

  delta = std::sqrt(delta);

  if (Math::are_float_equal(damage, 0.)) {
    delta_max_excess = 0;
  } else if (Math::are_float_equal(damage, 1.)) {
    delta_max_excess = 0;
  } else {
    delta_max_excess = delta / delta_max;
  }
  return delta_max_excess;
}

} // namespace akantu
/* ------------------------------------------------------------------- */
#endif // AKANTU_MATERIAL_COHESIVE_LINEAR_INLINE_IMPL_HH_
