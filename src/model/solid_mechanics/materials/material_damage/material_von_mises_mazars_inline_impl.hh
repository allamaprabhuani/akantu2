/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_von_mises_mazars.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
template <Int dim, template <UInt> class Parent>
inline void MaterialVonMisesMazars<dim, Parent>::computeStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam,
    Real & Ehat) {
  Matrix<Real> epsilon(3, 3);
  epsilon.zero();

  epsilon.block<dim, dim>(0, 0) = Material::gradUToEpsilon(grad_u);

  Vector<Real, 3> Fdiag(3);
  epsilon.eig(Fdiag);

  Ehat = 0.;
  for (UInt i = 0; i < 3; ++i) {
    Real epsilon_p = std::max(Real(0.), Fdiag(i));
    Ehat += epsilon_p * epsilon_p;
  }
  Ehat = std::sqrt(Ehat);

  // MaterialElastic<dim>::computeStressOnQuad(grad_u, sigma);

  if (damage_in_compute_stress) {
    computeDamageOnQuad(Ehat, sigma, Fdiag, dam);
  }

  if (not this->is_non_local) {
    computeDamageAndStressOnQuad(grad_u, sigma, dam, Ehat);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <UInt> class Parent>
inline void MaterialVonMisesMazars<dim, Parent>::computeDamageAndStressOnQuad(
    const Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam,
    Real & Ehat) {
  if (!damage_in_compute_stress) {
    auto && Fdiag = Vector<Real, 3>::Zero();

    auto && epsilon = Matrix<Real, 3, 3>::Zero();
    epsilon.block(0, 0, dim, dim) = Material::gradUToEpsilon<dim>(grad_u);
    epsilon.eig(Fdiag);

    computeDamageOnQuad(Ehat, sigma, Fdiag, dam);
  }

  sigma *= 1 - dam;
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <UInt> class Parent>
inline void MaterialVonMisesMazars<dim, Parent>::computeDamageOnQuad(
    const Real & epsilon_equ,
    __attribute__((unused)) const Matrix<Real> & sigma,
    const Vector<Real> & epsilon_princ, Real & dam) {
  Real Fs = epsilon_equ - K0;
  if (Fs > 0.) {
    Real dam_t;
    Real dam_c;
    dam_t =
        1 - K0 * (1 - At) / epsilon_equ - At * (exp(-Bt * (epsilon_equ - K0)));
    dam_c =
        1 - K0 * (1 - Ac) / epsilon_equ - Ac * (exp(-Bc * (epsilon_equ - K0)));

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

    alpha_t /= epsilon_equ * epsilon_equ;
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
