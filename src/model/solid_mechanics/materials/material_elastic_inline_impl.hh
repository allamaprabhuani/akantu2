/**
 * @file   material_elastic_inline_impl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Implementation of the inline functions of the material elastic
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
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__
// #define __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialElastic<dim>::computeStressOnQuad(Args && args) const {
  auto && sigma = tuple::get<"sigma"_h>(args);
  auto && grad_u = tuple::get<"grad_u"_h>(args);
  Real sigma_th = 0.;

  static_if(tuple::has_t<"sigma_th"_h, Args>()).then([&sigma_th](auto && args) {
    sigma_th = tuple::get<"sigma_th"_h>(args);
  })(std::forward<Args>(args));

  Real trace = grad_u.trace(); // trace = (\nabla u)_{kk}

  // \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  // u_{ij} + \nabla u_{ji})
  sigma = mu * Material::gradUToEpsilon<dim>(grad_u) +
          (lambda * trace + sigma_th) * Matrix<Real, dim, dim>::Identity();
}

/* -------------------------------------------------------------------------- */
template <>
template <typename Args>
inline void MaterialElastic<1>::computeStressOnQuad(Args && args) const {
  auto && sigma = tuple::get<"sigma"_h>(args);
  auto && grad_u = tuple::get<"grad_u"_h>(args);
  Real sigma_th = 0.;
  static_if(tuple::has_t<"sigma_th"_h, Args>()).then([&sigma_th](auto && args) {
    sigma_th = tuple::get<"sigma_th"_h>(args);
  })(std::forward<Args>(args));

  sigma(0, 0) = this->E * grad_u(0, 0) + sigma_th;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void
MaterialElastic<dim>::computeTangentModuliOnQuad(Args && args) const {
  auto && tangent = tuple::get<"tangent_moduli"_h>(args);
  constexpr auto n = Material::getTangentStiffnessVoigtSize(dim);

  // Real Ep = E/((1+nu)*(1-2*nu));
  auto Miiii = lambda + 2 * mu;
  auto Miijj = lambda;
  auto Mijij = mu;

  if (dim == 1) {
    tangent(0, 0) = this->E;
  } else {
    tangent(0, 0) = Miiii;
  }

  // test of dimension should by optimized out by the compiler due to the
  // template
  if (dim >= 2) {
    tangent(1, 1) = Miiii;
    tangent(0, 1) = Miijj;
    tangent(1, 0) = Miijj;

    tangent(n - 1, n - 1) = Mijij;
  }

  if (dim == 3) {
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
// template <>
// template <typename Derived>
// inline void MaterialElastic<1>::computeTangentModuliOnQuad(
//     Eigen::MatrixBase<Derived> & tangent) const {
//   tangent(0, 0) = E;
// }

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialElastic<dim>::computePotentialEnergyOnQuad(Args && args,
                                                               Real & epot) {
  epot =
      .5 * tuple::get<"sigma"_h>(args).doubleDot(tuple::get<"grad_u"_h>(args));
}

} // namespace akantu

// #endif /* __AKANTU_MATERIAL_ELASTIC_INLINE_IMPL_CC__ */
