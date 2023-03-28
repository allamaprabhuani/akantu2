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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialElastic<dim>::computeStressOnQuad(Args && args) const {
  auto && sigma = args["sigma"_n];
  auto && grad_u = args["grad_u"_n];
  Real sigma_th = 0.;

  if constexpr (named_tuple_t<Args>::has("sigma_th"_n)) {
    sigma_th = args["sigma_th"_n];
  }

  Real trace = grad_u.trace(); // trace = (\nabla u)_{kk}

  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})
  sigma = 2. * mu * Material::gradUToEpsilon<dim>(grad_u) +
          (lambda * trace + sigma_th) * Matrix<Real, dim, dim>::Identity();
}

/* -------------------------------------------------------------------------- */
template <>
template <typename Args>
inline void MaterialElastic<1>::computeStressOnQuad(Args && args) const {
  auto && sigma = args["sigma"_n];
  auto && grad_u = args["grad_u"_n];
  Real sigma_th = 0.;

  if constexpr (std::decay_t<Args>::has("sigma_th"_n)) {
    sigma_th = args["sigma_th"_n];
  }

  sigma(0, 0) = this->E * grad_u(0, 0) + sigma_th;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void
MaterialElastic<dim>::computeTangentModuliOnQuad(Args && args) const {
  auto && tangent = args["tangent_moduli"_n];
  tangent.zero();

  constexpr auto n = Material::getTangentStiffnessVoigtSize(dim);

  // Real Ep = E/((1+nu)*(1-2*nu));

  if constexpr (dim == 1) {
    tangent(0, 0) = this->E;
    return;
  }

  auto Miiii = lambda + 2 * mu;
  [[maybe_unused]] auto Miijj = lambda;
  [[maybe_unused]] auto Mijij = mu;

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

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElastic<dim>::computePotentialEnergyOnQuad(Args && args,
                                                               Real & epot) {
  epot = .5 * args["sigma"_n].doubleDot(args["grad_u"_n]);
}

} // namespace akantu

// #endif /* __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__ */
