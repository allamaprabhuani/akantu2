/**
 * @file   material_cohesive_extrinsic_bilinear_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Shenghan Zhang <shenghan.zhang@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Mar 5 2021
 * @date last modification: Fri Mar 5 2021
 *
 * @brief  Inline functions of the MaterialCohesiveExtrinsicBilinear
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
#include "material_cohesive_bilinear.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {


/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialCohesiveExtrinsicBilinear<dim>::computeTractionOnQuad(
    Vector<Real> & traction, Vector<Real> & opening,
    const Vector<Real> & normal, Real & delta_max, const Real & delta_c,
    const Real & delta_h, const Vector<Real> & insertion_stress,
    const Real & sigma_c, Vector<Real> & normal_opening,
    Vector<Real> & tangential_opening,
    Real & normal_opening_norm, Real & tangential_opening_norm, Real & damage,
    bool & penetration, Vector<Real> & contact_traction,
    Vector<Real> & contact_opening) {

  /// compute normal and tangential opening vectors
  normal_opening_norm = opening.dot(normal);
  normal_opening = (normal);
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
  
  penetration = normal_opening_norm / delta_c < -Math::getTolerance();
  if (not this->contact_after_breaking and
      Math::are_float_equal(damage, 1.)) {
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
  } else if (Math::are_float_equal(damage, 0)) {
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
    
    Real traction_norm;
    if (delta_max < delta_h) {
      Real sigma_h = h * delta_h;
      traction_norm = sigma_h + (delta_h - delta_max) / delta_h * (sigma_c - sigma_h);
    } else {
      Real sigma_h = h * delta_h;
      traction_norm = ((delta_c - delta_max) / (delta_c - delta_h)) * sigma_h;
    }
    
    traction *= traction_norm / delta_max;
  }
}


/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#endif //AKANTU_MATERIAL_COHESIVE_EXTRINSIC_BILINEAR_INLINE_IMPL_HH_
