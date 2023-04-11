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
#include "material_damage_iterative.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void
MaterialDamageIterative<spatial_dimension>::computeDamageAndStressOnQuad(
    Matrix<Real> & sigma, Real & dam) {
  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
UInt MaterialDamageIterative<spatial_dimension>::updateDamage(
    UInt quad_index, const Real /*eq_stress*/, ElementType el_type,
    GhostType ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  Array<Real> & dam = this->damage(el_type, ghost_type);
  Real & dam_on_quad = dam(quad_index);

  /// check if damage occurs
  if (equivalent_stress(el_type, ghost_type)(quad_index) >=
      (1 - dam_tolerance) * norm_max_equivalent_stress) {
    if (dam_on_quad < dam_threshold)
      dam_on_quad += prescribed_dam;
    else
      dam_on_quad = max_damage;
    return 1;
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline UInt MaterialDamageIterative<spatial_dimension>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {

  if (tag == SynchronizationTag::_user_2) {
    return sizeof(Real) * this->getModel().getNbIntegrationPoints(elements);
  }

  return MaterialDamage<spatial_dimension>::getNbData(elements, tag);
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void MaterialDamageIterative<spatial_dimension>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_user_2) {
    DataAccessor<Element>::packElementalDataHelper(
        this->damage, buffer, elements, true, this->damage.getFEEngine());
  }

  return MaterialDamage<spatial_dimension>::packData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
inline void MaterialDamageIterative<spatial_dimension>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_user_2) {
    DataAccessor<Element>::unpackElementalDataHelper(
        this->damage, buffer, elements, true, this->damage.getFEEngine());
  }
  return MaterialDamage<spatial_dimension>::unpackData(buffer, elements, tag);
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
