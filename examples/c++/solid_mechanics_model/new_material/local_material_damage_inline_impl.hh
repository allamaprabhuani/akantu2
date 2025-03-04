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
#ifndef AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH_
#define AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
inline void LocalMaterialDamage::computeStressOnQuad(
    Eigen::MatrixBase<D1> &grad_u, Eigen::MatrixBase<D2> &sigma, Real &dam) {

  Real trace = grad_u.trace();

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  /// u_{ij} + \nabla u_{ji})
  for (Int i = 0; i < spatial_dimension; ++i) {
    for (Int j = 0; j < spatial_dimension; ++j) {
      sigma(i, j) =
          (i == j) * lambda * trace + mu * (grad_u(i, j) + grad_u(j, i));
    }
  }

  Real Y = 0;
  for (Int i = 0; i < spatial_dimension; ++i) {
    for (Int j = 0; j < spatial_dimension; ++j) {
      Y += sigma(i, j) * grad_u(i, j);
    }
  }
  Y *= 0.5;

  Real Fd = Y - Yd - Sd * dam;

  if (Fd > 0)
    dam = (Y - Yd) / Sd;
  dam = std::min(dam, 1.);

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <typename D1, typename D2>
inline void LocalMaterialDamage::computePotentialEnergyOnQuad(
    Eigen::MatrixBase<D1> &grad_u, Eigen::MatrixBase<D2> &sigma, Real &epot) {
  epot = 0.;
  for (Int i = 0, t = 0; i < spatial_dimension; ++i)
    for (Int j = 0; j < spatial_dimension; ++j, ++t)
      epot += sigma(i, j) * (grad_u(i, j) - (i == j));
  epot *= .5;
}

/* -------------------------------------------------------------------------- */
inline Real LocalMaterialDamage::getCelerity(__attribute__((unused))
                                             const Element &element) const {
  // Here the fastest celerity is the push wave speed
  return (std::sqrt((2 * mu + lambda) / rho));
}

} // namespace akantu

#endif /* AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH_ */
