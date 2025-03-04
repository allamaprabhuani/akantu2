/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_iterators.hh"
#include "material_anisotropic_damage.hh" // NOLINT(pp_including_mainfile_in_preamble)
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_TMPL_HH_
#define AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_TMPL_HH_

namespace akantu {
struct EmptyIteratorContainer {
  using size_type = std::size_t;

  struct iterator {
    auto & operator++() { return *this; }
    Real operator*() { return 0; }
    bool operator!=(const iterator & /*unused*/) const { return true; }
    bool operator==(const iterator & /*unused*/) const { return false; }
  };

  auto begin() const // NOLINT(readability-convert-member-functions-to-static)
  {
    return iterator();
  }

  auto end() const // NOLINT(readability-convert-member-functions-to-static)
  {
    return iterator();
  }
};
} // namespace akantu

namespace std {
template <> struct iterator_traits<::akantu::EmptyIteratorContainer::iterator> {
  using iterator_category = forward_iterator_tag;
  using value_type = akantu::Real;
  using difference_type = std::ptrdiff_t;
  using pointer = akantu::Real *;
  using reference = akantu::Real &;
};
} // namespace std

namespace akantu {
namespace {
  template <Int dim, class D, class Op>
  void tensorPlus_(const Eigen::MatrixBase<D> & A, Op && oper) {
    Vector<Real, dim> A_eigs;
    A.eig(A_eigs);

    for (auto & ap : A_eigs) {
      oper(ap);
    }
  }

  template <Int dim, class D> auto tensorPlus2(const Eigen::MatrixBase<D> & A) {
    Real square = 0;
    tensorPlus_<dim>(A, [&](Real eig) {
      eig = std::max(eig, 0.);
      square += eig * eig;
    });

    return square;
  }

  template <Int dim, class D>
  auto tensorPlusTrace(const Eigen::MatrixBase<D> & A) {
    Real trace_plus = 0.;
    Real trace_minus = 0.;
    tensorPlus_<dim>(A, [&](Real eig) {
      trace_plus += std::max(eig, 0.);
      trace_minus += std::min(eig, 0.);
    });

    return std::pair(trace_plus, trace_minus);
  }

  template <Int dim, class Op, class D1, class D2>
  auto tensorPlusOp(const Eigen::MatrixBase<D1> & A,
                    Eigen::MatrixBase<D2> & A_directions, Op && oper) {
    Vector<Real, dim> A_eigs;
    Matrix<Real, dim, dim> A_diag;
    A.eig(A_eigs, A_directions);

    for (auto && [i, eig] : enumerate(A_eigs)) {
      A_diag(i, i) = oper(std::max(eig, 0.), i);
    }

    Matrix<Real, dim, dim> res =
        A_directions * A_diag * A_directions.transpose();
    return res;
  }

  template <Int dim, class D1, class D2, class Op>
  auto tensorPlus(const Eigen::MatrixBase<D1> & A,
                  Eigen::MatrixBase<D2> & A_directions) {
    return tensorPlusOp<dim>(A, A_directions, [](Real x, Real) { return x; });
  }

  template <Int dim, class D, class Op>
  auto tensorPlusOp(const Eigen::MatrixBase<D> & A, Op && oper) {
    Matrix<Real, dim, dim> A_directions;
    A_directions.zero();
    return tensorPlusOp<dim>(A, A_directions, std::forward<Op>(oper));
  }

  template <Int dim, class D> auto tensorPlus(const Eigen::MatrixBase<D> & A) {
    return tensorPlusOp<dim>(A, [](Real x, Real) { return x; });
  }

  template <Int dim, class D> auto tensorSqrt(const Eigen::MatrixBase<D> & A) {
    return tensorPlusOp<dim>(A, [](Real x, Int) { return std::sqrt(x); });
  }

} // namespace

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EquivalentStrain,
          template <Int> class DamageThreshold, template <Int> class Parent>
MaterialAnisotropicDamage<dim, EquivalentStrain, DamageThreshold, Parent>::
    MaterialAnisotropicDamage(SolidMechanicsModel & model, const ID & id)
    : Parent<dim>(model, id), equivalent_strain_function(*this),
      damage_threshold_function(*this),
      damage(this->registerInternal("damage", dim * dim)),
      elastic_stress(this->registerInternal("elastic_stress", dim * dim)),
      equivalent_strain(this->registerInternal("equivalent_strain", 1)),
      trace_damage(this->registerInternal("trace_damage", 1)) {
  this->registerParam("Dc", Dc, _pat_parsable, "Critical damage");
  this->trace_damage.initializeHistory();
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EquivalentStrain,
          template <Int> class DamageThreshold, template <Int> class Parent>
template <class D1, class D2, class D3>
void MaterialAnisotropicDamage<dim, EquivalentStrain, DamageThreshold, Parent>::
    damageStress(Eigen::MatrixBase<D1> & sigma,
                 const Eigen::MatrixBase<D2> & sigma_el,
                 const Eigen::MatrixBase<D3> & D, Real TrD) {
  // σ_(n + 1) = (1 − D_(n + 1))^(1/2) σ~_(n + 1) (1 − D_(n + 1))^(1 / 2)
  //         - ((1 − D_(n + 1)) : σ~_(n + 1))/ (3 - Tr(D_(n+1))) (1 − D_(n + 1))
  //         + 1/3 (1 - Tr(D_(n+1)) <Tr(σ~_(n + 1))>_+ + <Tr(σ~_(n + 1))>_-) I

  Matrix<Real, dim, dim> one_D = Matrix<Real, dim, dim>::Identity() - D;
  auto sqrt_one_D = tensorSqrt<dim>(one_D);

  auto && [Tr_sigma_plus, Tr_sigma_minus] = tensorPlusTrace<dim>(sigma_el);
  auto I = Matrix<Real, dim, dim>::Identity();

  sigma = sqrt_one_D * sigma_el * sqrt_one_D -
          (one_D.doubleDot(sigma_el) / (dim - TrD) * one_D) +
          1. / dim * ((1 - TrD) * Tr_sigma_plus - Tr_sigma_minus) * I;
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EquivalentStrain,
          template <Int> class DamageThreshold, template <Int> class Parent>
template <class Args>
void MaterialAnisotropicDamage<dim, EquivalentStrain, DamageThreshold,
                               Parent>::computeStressOnQuad(Args && args) {
  auto & sigma = args["sigma"_n];
  auto & grad_u = args["grad_u"_n];
  auto & sigma_el = args["sigma_el"_n];
  auto & epsilon_hat = args["epsilon_hat"_n];
  auto & D = args["damage"_n];
  auto & TrD_n_1 = args["TrD_n_1"_n];
  auto & TrD = args["TrD"_n];
  auto & equivalent_strain_data = args["equivalent_strain_data"_n];
  auto & damage_threshold_data = args["damage_threshold_data"_n];

  Matrix<Real, dim, dim> Dtmp;
  Real TrD_n_1_tmp{0.};
  Matrix<Real, dim, dim> epsilon;

  // yes you read properly this is a label for a goto
  auto computeDamage = [&]() {
    MaterialElastic<dim>::computeStressOnQuad(args);

    epsilon = this->template gradUToEpsilon<dim>(grad_u);

    // evaluate the damage criteria
    epsilon_hat = equivalent_strain_function(epsilon, equivalent_strain_data);

    // evolve the damage if needed
    auto K_TrD = damage_threshold_function.K(TrD, damage_threshold_data);

    auto f = epsilon_hat - K_TrD;

    // if test function > 0 evolve the damage
    if (f > 0) {
      TrD_n_1_tmp =
          damage_threshold_function.K_inv(epsilon_hat, damage_threshold_data);

      auto epsilon_plus = tensorPlus<dim>(epsilon);
      auto delta_lambda = (TrD_n_1_tmp - TrD) / (epsilon_hat * epsilon_hat);

      Dtmp = D + delta_lambda * epsilon_plus;
      return true;
    }
    return false;
  };

  // compute a temporary version of the new damage
  auto is_damage_updated = computeDamage();

  if (is_damage_updated) {
    /// Check and correct for broken case
    if (Dtmp.trace() > Dc) {
      if (epsilon.trace() > 0) { // tensile loading
        auto kpa = this->kpa;
        auto lambda = this->lambda;

        // change kappa to Kappa_broken = (1-Dc) Kappa
        kpa = (1 - Dc) * kpa;
        this->E = 9 * kpa * (kpa - lambda) / (3 * kpa - lambda);
        this->nu = lambda / (3 * kpa - lambda);
        this->updateInternalParameters();

        computeDamage();
      } else if (std::abs(epsilon.trace()) < 1e-10) { // deviatoric case
        Matrix<Real, dim, dim> n;
        std::vector<Int> ns;
        tensorPlusOp<dim>(Dtmp, n, [&](Real x, Int i) {
          if (x > this->Dc) {
            ns.push_back(i);
            return this->Dc;
          }

          return x;
        });
      }
    }

    TrD_n_1 = TrD_n_1_tmp;
    D = Dtmp;
  } else {
    TrD_n_1 = TrD;
  }

  // apply the damage to the stress
  damageStress(sigma, sigma_el, D, TrD_n_1);
}

/* -------------------------------------------------------------------------- */
template <Int dim, template <Int> class EquivalentStrain,
          template <Int> class DamageThreshold, template <Int> class Parent>
void MaterialAnisotropicDamage<dim, EquivalentStrain, DamageThreshold,
                               Parent>::computeStress(ElementType type,
                                                      GhostType ghost_type) {
  for (auto && args : getArguments(type, ghost_type)) {
    computeStressOnQuad(args);
  }
}

/* -------------------------------------------------------------------------- */
/* EquivalentStrain functions                                                 */
/* -------------------------------------------------------------------------- */
template <Int dim>
class EquivalentStrainMazars : public EmptyIteratorContainer {
public:
  EquivalentStrainMazars(Material & /*mat*/) {}

  template <class D, class... Other>
  Real operator()(const Eigen::MatrixBase<D> & epsilon, Other &&... /*other*/) {
    auto && [epsilon_hat, _] = tensorPlusTrace<dim>(epsilon);
    return std::sqrt(epsilon_hat);
  }
};

template <Int dim>
class EquivalentStrainMazarsDruckerPrager : public EquivalentStrainMazars<dim> {
public:
  EquivalentStrainMazarsDruckerPrager(Material & mat)
      : EquivalentStrainMazars<dim>(mat) {
    mat.registerParam("k", k, _pat_parsable, "k");
  }

  template <class D, class... Other>
  Real operator()(const Eigen::MatrixBase<D> & epsilon, Other &&... /*other*/) {
    Real epsilon_hat = EquivalentStrainMazars<dim>::operator()(epsilon);
    epsilon_hat += k * epsilon.trace();
    return epsilon_hat;
  }

protected:
  Real k{0.};
};

/* -------------------------------------------------------------------------- */
/* DamageThreshold functions                                                  */
/* -------------------------------------------------------------------------- */
template <Int dim> class DamageThresholdLinear : public EmptyIteratorContainer {
public:
  DamageThresholdLinear(Material & mat) : mat(mat) {
    mat.registerParam("A", A, _pat_parsable, "A");
    mat.registerParam("K0", K0, _pat_parsable, "K0");
  }

  template <class... Other> Real K(Real x, Other &&... /*other*/) {
    return 1. / A * x + K0;
  }

  template <class... Other> Real K_inv(Real x, Other &&... /*other*/) {
    return A * (x - K0);
  }

private:
  Material & mat;
  Real A{0.};
  Real K0{0.};
};

template <Int dim> class DamageThresholdTan : public EmptyIteratorContainer {
public:
  DamageThresholdTan(Material & mat) : mat(mat) {
    mat.registerParam("a", a, _pat_parsable, "a");
    mat.registerParam("A", A, _pat_parsable, "A");
    mat.registerParam("K0", K0, _pat_parsable, "K0");
  }

  template <class... Other> Real K(Real x, Other &&... /*other*/) {
    return a * std::tan(std::atan2(x, a) - std::atan2(K0, a));
  }

  template <class... Other> Real K_inv(Real x, Other &&... /*other*/) {
    return a * A * (std::atan2(x, a) - std::atan2(K0, a));
  }

private:
  Material & mat;
  Real a{2.93e-4};
  Real A{5e3};
  Real K0{5e-5};
};

} // namespace akantu

#endif /* AKANTU_MATERIAL_ANISOTROPIC_DAMAGE_TMPL_HH_ */
