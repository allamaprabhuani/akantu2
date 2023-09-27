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
#include "local_material_damage.hh" // NOLINT
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__
#define AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class D1, class D2>
inline void LocalMaterialDamage::computeStressOnQuad(
    Eigen::MatrixBase<D1> & grad_u, Eigen::MatrixBase<D2> & sigma, Real & dam) {

  Real trace = grad_u.trace();

  /// \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  /// u_{ij} + \nabla u_{ji})
  auto && epsilon = (grad_u + grad_u.transpose()) / 2.;

  sigma = Matrix<Real>::Identity(spatial_dimension, spatial_dimension) * trace *
              lambda +
          mu * epsilon;

  Real Y = sigma.doubleDot(epsilon) / 2.;
  Real Fd = Y - Yd - Sd * dam;

  if (Fd > 0) {
    dam = (Y - Yd) / Sd;
  }
  dam = std::min(dam, 1.);

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <class D1, class D2>
inline void LocalMaterialDamage::computePotentialEnergyOnQuad(
    Eigen::MatrixBase<D1> & grad_u, Eigen::MatrixBase<D2> & sigma,
    Real & epot) {
  epot = sigma.doubleDot(grad_u) / 2.;
}

/* -------------------------------------------------------------------------- */
inline Real
LocalMaterialDamage::getCelerity(const Element & /*element*/) const {
  return (std::sqrt(E / rho));
}

} // namespace akantu

#endif /* AKANTU_LOCAL_MATERIAL_DAMAGE_INLINE_IMPL_HH__ */
