/**
 * @file   material_phasefield_anisotropic_inline_impl.cc
 *
 * @author Shad Durussel <shad.durussel@epfl.ch>
 *
 * @date creation: Mon Mar 27 2023
 * @date last modification: Mon Mar 27 2023
 *
 * @brief  Implementation of the inline functions of the material phasefield
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

#include "material_phasefield_anisotropic.hh"
#include <algorithm>

#ifndef AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_INLINE_IMPL_HH_
#define AKANTU_MATERIAL_PHASEFIELD_ANISOTROPIC_INLINE_IMPL_HH_
/* -------------------------------------------------------------------------- */
namespace akantu {
template <Int dim>
template <class Args>
inline void
MaterialPhaseFieldAnisotropic<dim>::computeStressOnQuad(Args && args) {
  // MaterialElastic<dim>::computeStressOnQuad(args);

  auto && dam = args["damage"_n];

  auto sigma = args["sigma"_n];
  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  Real trace = strain.trace();
  Real trace_plus = std::max(Real(0.), trace);
  Real trace_minus = std::min(Real(0.), trace);

  Real sigma_th_plus = std::max(Real(0.), args["sigma_th"_n]);
  Real sigma_th_minus = std::min(Real(0.), args["sigma_th"_n]);

  Matrix<Real> strain_dev(dev_dim, dev_dim);
  Matrix<Real> strain_tmp = Matrix<Real>::Zero(dev_dim, dev_dim);
  strain_tmp.topLeftCorner(dim, dim) = strain;

  strain_dev = strain_tmp - trace / Real(dev_dim) * Matrix<Real>::Identity(dev_dim, dev_dim);

  Real kappa = this->lambda + 2. / Real(dev_dim) * this->mu;

  Real g_d = (1 - dam) * (1 - dam) + eta;

  auto sigma_plus = (kappa * trace_plus + sigma_th_plus) *
                        Matrix<Real>::Identity(dev_dim, dev_dim) +
                    2. * this->mu * strain_dev;
  auto sigma_minus = (kappa * trace_minus + sigma_th_minus) *
                     Matrix<Real>::Identity(dev_dim , dev_dim);

  sigma = g_d * sigma_plus.topLeftCorner(dim, dim) + sigma_minus.topLeftCorner(dim, dim);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
void MaterialPhaseFieldAnisotropic<dim>::computeTangentModuliOnQuad(
    Args && args) {

  auto && dam = args["damage"_n];

  auto strain = Material::gradUToEpsilon<dim>(args["grad_u"_n]);

  Real trace = strain.trace();

  Real g_d = (1 - dam) * (1 - dam) + eta;
  Real g_d_hyd = trace > 0. ? g_d : 1.;

  auto && tangent = args["tangent_moduli"_n];
  tangent.zero();

  constexpr auto n = Material::getTangentStiffnessVoigtSize(dim);

  // Real Ep = E/((1+nu)*(1-2*nu));

  if constexpr (dim == 1) {
    tangent(0, 0) = g_d_hyd * this->E;
    return;
  }

  Real kappa = this->lambda + 2. / Real(dev_dim) * this->mu;

  auto Miiii = g_d_hyd * kappa + g_d * 2. * this->mu * (1. - 1. / Real(dev_dim));
  [[maybe_unused]] auto Miijj = g_d_hyd * kappa - g_d * 2. * this->mu / Real(dev_dim);
  [[maybe_unused]] auto Mijij = g_d * this->mu;

  tangent(0, 0) = Miiii;

  // test of dimension should by optimized out by the compiler due to the
  // template
  if constexpr (dim >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if constexpr (dim == 3) {
    tangent(2, 2) = Miiii;
    tangent(0, 2) = Miijj;
    tangent(1, 2) = Miijj;
    tangent(2, 0) = Miijj;
    tangent(2, 1) = Miijj;

    tangent(3, 3) = Mijij;
    tangent(4, 4) = Mijij;
  }
}

} // namespace akantu
#endif
