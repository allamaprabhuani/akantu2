/**
 * @file   material_marigo_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implementation of the inline functions of the material marigo
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_marigo.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialMarigo<dim>::computeStressOnQuad(Args && arguments) {
  auto && sigma = tuple::get<"sigma"_h>(arguments);
  auto && grad_u = tuple::get<"grad_u"_h>(arguments);
  auto && dam = tuple::get<"damage"_h>(arguments);
  auto && Y = tuple::get<"Y"_h>(arguments);

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
  auto && sigma = tuple::get<"sigma"_h>(arguments);
  auto && dam = tuple::get<"damage"_h>(arguments);
  auto && Y = tuple::get<"Y"_h>(arguments);
  auto && Yd = tuple::get<"Yd"_h>(arguments);

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

  Int size = 0;
  if (tag == SynchronizationTag::_smm_init_mat) {
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

  if (tag == SynchronizationTag::_smm_init_mat) {
    this->packElementDataHelper(Yd, buffer, elements);
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

  if (tag == SynchronizationTag::_smm_init_mat) {
    this->unpackElementDataHelper(Yd, buffer, elements);
  }

  MaterialDamage<dim>::unpackData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

//#endif /* __AKANTU_MATERIAL_MARIGO_INLINE_IMPL_CC__ */
