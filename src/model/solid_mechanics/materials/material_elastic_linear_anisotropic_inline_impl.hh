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
  auto && sigma = args["sigma"_n];
  auto && grad_u = args["grad_u"_n];

  auto voigt_strain = strainToVoigt<dim>(gradUToEpsilon<dim>(grad_u));
  auto voigt_stress = this->C * voigt_strain;
  voigtToStress<dim>(voigt_stress, sigma);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElasticLinearAnisotropic<dim>::computeTangentModuliOnQuad(
    Args && args) const {
  args["tangent_moduli"_n] = this->C;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElasticLinearAnisotropic<dim>::computePotentialEnergyOnQuad(
    Args && args, Real & epot) {

  AKANTU_DEBUG_ASSERT(this->symmetric,
                      "The elastic constants matrix is not symmetric,"
                      "energy is not path independent.");

  epot = args["sigma"_n].doubleDot(args["grad_u"_n]) / 2.;
}

} // namespace akantu

#endif /* AKANTU_MATERIAL_ELASTIC_LINEAR_ANISOTROPIC_INLINE_IMPL_HH_ */
