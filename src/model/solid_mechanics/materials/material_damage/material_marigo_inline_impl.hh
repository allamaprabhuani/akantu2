/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_marigo.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialMarigo<dim>::computeStressOnQuad(Args && arguments) {
  auto && sigma = arguments["sigma"_n];
  auto && grad_u = arguments["grad_u"_n];
  auto && dam = arguments["damage"_n];
  auto && Y = arguments["Y"_n];

  MaterialElastic<dim>::computeStressOnQuad(arguments);

  Y = sigma.doubleDot(Material::gradUToEpsilon<dim>(grad_u)) / 2.;

  if (this->damage_in_y) {
    Y *= (1 - dam);
  }

  if (this->yc_limit) {
    Y = std::min(Y, this->Yc);
  }

  if (not this->is_non_local) {
    computeDamageAndStressOnQuad(arguments);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void
MaterialMarigo<dim>::computeDamageAndStressOnQuad(Args && arguments) {
  auto && sigma = arguments["sigma"_n];
  auto && dam = arguments["damage"_n];
  auto && Y = arguments["Y"_n];
  auto && Yd = arguments["Yd"_n];

  Real Fd = Y - Yd - Sd * dam;

  if (Fd > 0) {
    dam = (Y - Yd) / Sd;
  }
  dam = std::min(dam, Real(1.));

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline auto MaterialMarigo<dim>::getNbData(const Array<Element> & elements,
                                           const SynchronizationTag & tag) const
    -> Int {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if (tag == SynchronizationTag::_clh_init_cl) {
    size += sizeof(Real) * this->getModel().getNbIntegrationPoints(elements);
  }

  size += MaterialDamage<dim>::getNbData(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
MaterialMarigo<dim>::packData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  if (tag == SynchronizationTag::_clh_init_cl) {
    this->packInternalFieldHelper(Yd, buffer, elements);
  }

  MaterialDamage<dim>::packData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void MaterialMarigo<dim>::unpackData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (tag == SynchronizationTag::_clh_init_cl) {
    this->unpackInternalFieldHelper(Yd, buffer, elements);
  }

  MaterialDamage<dim>::unpackData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

//#endif /* __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__ */
