/**
 * @file   material_elastic_linear_anisotropic_inline_impl.hh
 *
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Feb 16 2018
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Implementation of the inline functions of the material elastic linear
 * anisotropic
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_elastic_linear_anisotropic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void
MaterialElasticLinearAnisotropic<dim>::computeStressOnQuad(Args && args) const {
  auto && sigma = tuple::get<"sigma"_h>(args);
  auto && grad_u = tuple::get<"grad_u"_h>(args);

  auto voigt_strain = strainToVoigt<dim>(gradUToEpsilon<dim>(grad_u));
  auto voigt_stress = this->C * voigt_strain;
  voigtToStress<dim>(voigt_stress, sigma);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElasticLinearAnisotropic<dim>::computeTangentModuliOnQuad(
    Args && args) const {
  tuple::get<"tangent_moduli"_h>(args) = this->C;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergyOnQuad(
    Args && args, Real & epot) {

  AKANTU_DEBUG_ASSERT(this->symmetric,
                      "The elastic constants matrix is not symmetric,"
                      "energy is not path independent.");

  epot =
      tuple::get<"sigma"_h>(args).doubleDot(tuple::get<"grad_u"_h>(args)) / 2.;
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_ */
