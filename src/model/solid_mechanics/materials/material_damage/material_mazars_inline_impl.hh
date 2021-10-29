/**
 * @file   material_mazars_inline_impl.hh
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Marion Estelle Chambart <mchambart@stucky.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 06 2011
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Implementation of the inline functions of the material damage
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
#include "material_mazars.hh"

namespace akantu {
/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void MaterialMazars<dim>::computeStressOnQuad(Args && arguments) {
  auto && grad_u = tuple::get<"grad_u"_h>(arguments);
  auto && Ehat = tuple::get<"Ehat"_h>(arguments);

  Matrix<Real, 3, 3> epsilon;
  epsilon.fill(0.);

  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      epsilon(i, j) = .5 * (grad_u(i, j) + grad_u(j, i));
    }
  }

  Vector<Real, 3> Fdiag;
  epsilon.eigh(Fdiag);

  Ehat = 0.;
  for (Int i = 0; i < 3; ++i) {
    Real epsilon_p = std::max(Real(0.), Fdiag(i));
    Ehat += epsilon_p * epsilon_p;
  }
  Ehat = sqrt(Ehat);

  MaterialElastic<dim>::computeStressOnQuad(arguments);

  if (damage_in_compute_stress) {
    computeDamageOnQuad(arguments, Fdiag);
  }

  if (not this->is_non_local) {
    computeDamageAndStressOnQuad(arguments);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args>
inline void
MaterialMazars<dim>::computeDamageAndStressOnQuad(Args && arguments) {
  auto && sigma = tuple::get<"sigma"_h>(arguments);
  auto && grad_u = tuple::get<"grad_u"_h>(arguments);
  auto && dam = tuple::get<"damage"_h>(arguments);

  if (!damage_in_compute_stress) {
    Vector<Real, 3> Fdiag;
    Fdiag.clear();

    Matrix<Real, 3, 3> epsilon;
    epsilon.clear();
    for (Int i = 0; i < dim; ++i) {
      for (Int j = 0; j < dim; ++j) {
        epsilon(i, j) = .5 * (grad_u(i, j) + grad_u(j, i));
      }
    }

    epsilon.eigh(Fdiag);

    computeDamageOnQuad(arguments, Fdiag);
  }

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <typename Args, typename Derived>
inline void MaterialMazars<dim>::computeDamageOnQuad(
    Args && arguments, const Eigen::MatrixBase<Derived> & epsilon_princ) {
  auto && dam = tuple::get<"damage"_h>(arguments);
  auto && Ehat = tuple::get<"Ehat"_h>(arguments);

  Real Fs = Ehat - K0;
  if (Fs > 0.) {
    Real dam_t;
    Real dam_c;
    dam_t =
        1 - K0 * (1 - At) / Ehat - At * (exp(-Bt * (Ehat - K0)));
    dam_c =
        1 - K0 * (1 - Ac) / Ehat - Ac * (exp(-Bc * (Ehat - K0)));

    Real Cdiag;
    Cdiag = this->E * (1 - this->nu) / ((1 + this->nu) * (1 - 2 * this->nu));

    Vector<Real> sigma_princ(3);
    sigma_princ(0) = Cdiag * epsilon_princ(0) +
                     this->lambda * (epsilon_princ(1) + epsilon_princ(2));
    sigma_princ(1) = Cdiag * epsilon_princ(1) +
                     this->lambda * (epsilon_princ(0) + epsilon_princ(2));
    sigma_princ(2) = Cdiag * epsilon_princ(2) +
                     this->lambda * (epsilon_princ(1) + epsilon_princ(0));

    Vector<Real> sigma_p(3);
    for (Int i = 0; i < 3; i++) {
      sigma_p(i) = std::max(Real(0.), sigma_princ(i));
    }
    // sigma_p *= 1. - dam;

    Real trace_p = this->nu / this->E * (sigma_p(0) + sigma_p(1) + sigma_p(2));

    Real alpha_t = 0;
    for (Int i = 0; i < 3; ++i) {
      Real epsilon_t = (1 + this->nu) / this->E * sigma_p(i) - trace_p;
      Real epsilon_p = std::max(Real(0.), epsilon_princ(i));
      alpha_t += epsilon_t * epsilon_p;
    }

    alpha_t /= Ehat * Ehat;
    alpha_t = std::min(alpha_t, Real(1.));

    Real alpha_c = 1. - alpha_t;

    alpha_t = std::pow(alpha_t, beta);
    alpha_c = std::pow(alpha_c, beta);

    Real damtemp;
    damtemp = alpha_t * dam_t + alpha_c * dam_c;

    dam = std::max(damtemp, dam);
    dam = std::min(dam, Real(1.));
  }
}

} // namespace akantu
