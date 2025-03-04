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
// #include "material_linear_isotropic_hardening.hh"
#include "material_mazars.hh" // NOLINT

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
MaterialMazars<dim, Parent>::MaterialMazars(SolidMechanicsModel & model,
                                            const ID & id)
    : parent_damage(model, id),
      K0(this->template registerInternal<Real, DefaultRandomInternalField>("K0",
                                                                           1)) {

  this->registerParam("K0", this->K0, _pat_parsable, "K0");
  this->registerParam("At", this->At, Real(0.8), _pat_parsable, "At");
  this->registerParam("Ac", this->Ac, Real(1.4), _pat_parsable, "Ac");
  this->registerParam("Bc", this->Bc, Real(1900.), _pat_parsable, "Bc");
  this->registerParam("Bt", this->Bt, Real(12000.), _pat_parsable, "Bt");
  this->registerParam("beta", this->beta, Real(1.06), _pat_parsable, "beta");
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
void MaterialMazars<dim, Parent>::computeStress(ElementType el_type,
                                                GhostType ghost_type) {
  auto && arguments = getArguments(el_type, ghost_type);
  for (auto && args : arguments) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
template <typename Args>
inline void MaterialMazars<dim, Parent>::computeStressOnQuad(Args && args) {
  Parent<dim>::computeStressOnQuad(std::forward<decltype(args)>(args));

  auto & grad_u = args["grad_u"_n];

  if constexpr (named_tuple_t<Args>::has("inelastic_strain"_n)) {
    grad_u -= args["inelastic_strain"_n];
  }

  Matrix<Real, 3, 3> epsilon = Matrix<Real, 3, 3>::Zero();

  epsilon.block<dim, dim>(0, 0) = Material::gradUToEpsilon<dim>(grad_u);

  Vector<Real, 3> Fdiag;
  epsilon.eig(Fdiag);

  auto & Ehat = args["Ehat"_n];

  Ehat = 0.;
  for (Int i = 0; i < 3; ++i) {
    Real epsilon_p = std::max(Real(0.), Fdiag(i));
    Ehat += epsilon_p * epsilon_p;
  }
  Ehat = std::sqrt(Ehat);

  if (damage_in_compute_stress) {
    computeDamageOnQuad(args, Fdiag);
  }

  if (not this->is_non_local) {
    computeDamageAndStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
template <typename Args>
inline void
MaterialMazars<dim, Parent>::computeDamageAndStressOnQuad(Args && args) {
  auto && grad_u = args["grad_u"_n];
  if (not damage_in_compute_stress) {
    Vector<Real, 3> Fdiag;
    Matrix<Real, 3, 3> epsilon = Matrix<Real, 3, 3>::Zero();

    epsilon.block<dim, dim>(0, 0) = Material::gradUToEpsilon<dim>(grad_u);

    epsilon.eig(Fdiag);
    computeDamageOnQuad(std::forward<Args>(args), Fdiag);
  }

  auto && sigma = args["sigma"_n];
  auto && dam = args["damage"_n];
  sigma *= 1 - dam;

  if constexpr (named_tuple_t<Args>::has("inelastic_strain"_n)) {
    grad_u += args["inelastic_strain"_n];
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class Parent>
template <typename Args, typename Derived>
inline void MaterialMazars<dim, Parent>::computeDamageOnQuad(
    Args && args, const Eigen::MatrixBase<Derived> & epsilon_princ) {
  auto && dam = args["damage"_n];
  auto && Ehat = args["Ehat"_n];
  auto && K0 = args["K0"_n];

  auto Fs = Ehat - K0;

  if (Fs <= 0.) {
    return;
  }

  auto dam_t = 1 - K0 * (1 - At) / Ehat - At * (exp(-Bt * (Ehat - K0)));
  auto dam_c = 1 - K0 * (1 - Ac) / Ehat - Ac * (exp(-Bc * (Ehat - K0)));

  auto Cdiag = this->E * (1 - this->nu) / ((1 + this->nu) * (1 - 2 * this->nu));

  Vector<Real, 3> sigma_princ;
  sigma_princ(0) = Cdiag * epsilon_princ(0) +
                   this->lambda * (epsilon_princ(1) + epsilon_princ(2));
  sigma_princ(1) = Cdiag * epsilon_princ(1) +
                   this->lambda * (epsilon_princ(0) + epsilon_princ(2));
  sigma_princ(2) = Cdiag * epsilon_princ(2) +
                   this->lambda * (epsilon_princ(1) + epsilon_princ(0));

  Vector<Real, 3> sigma_p;
  for (Int i = 0; i < 3; i++) {
    sigma_p(i) = std::max(Real(0.), sigma_princ(i));
  }
  // sigma_p *= 1. - dam;

  auto trace_p = this->nu / this->E * (sigma_p(0) + sigma_p(1) + sigma_p(2));

  Real alpha_t = 0;
  for (Int i = 0; i < 3; ++i) {
    auto epsilon_t = (1 + this->nu) / this->E * sigma_p(i) - trace_p;
    auto epsilon_p = std::max(Real(0.), epsilon_princ(i));
    alpha_t += epsilon_t * epsilon_p;
  }

  alpha_t /= Ehat * Ehat;
  alpha_t = std::min(alpha_t, Real(1.));

  auto alpha_c = 1. - alpha_t;

  alpha_t = std::pow(alpha_t, beta);
  alpha_c = std::pow(alpha_c, beta);

  auto damtemp = alpha_t * dam_t + alpha_c * dam_c;

  dam = std::max(damtemp, dam);
  dam = std::min(dam, Real(1.));
}

} // namespace akantu
