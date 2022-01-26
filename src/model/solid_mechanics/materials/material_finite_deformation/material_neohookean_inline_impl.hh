/**
 * @file   material_neohookean_inline_impl.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Mon Apr 08 2013
 * @date last modification: Thu Feb 20 2020
 *
 * @brief  Implementation of the inline functions of the material elastic
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_neohookean.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
#include <utility>
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
template <Int dim>
inline void MaterialNeohookean<dim>::computeDeltaStressOnQuad(
    __attribute__((unused)) const Matrix<Real> & grad_u,
    __attribute__((unused)) const Matrix<Real> & grad_delta_u,
    __attribute__((unused)) Matrix<Real> & delta_S) {}

//! computes the second piola kirchhoff stress, called S
template <Int dim>
template <class Args>
inline void MaterialNeohookean<dim>::computeStressOnQuad(Args && args) {
  // Neo hookean book
  auto && F = Material::gradUToF<dim>(tuple::get<"grad_u"_h>(args));
  auto && C = Material::rightCauchy<dim>(F);
  // the term  sqrt(C33) corresponds to the off plane strain (2D plane stress)
  auto J = F.determinant() * std::sqrt(tuple::get<"C33"_h>(args));

  tuple::get<"sigma"_h>(args) = Matrix<Real, dim, dim>::Identity() * mu +
                                (lambda * log(J) - mu) * C.inverse();
}

/* -------------------------------------------------------------------------- */
class C33_NR : public Math::NewtonRaphsonFunctor<Real> {
public:
  C33_NR(std::string name, const Real & lambda, const Real & mu,
         const Matrix<Real> & C)
      : NewtonRaphsonFunctor(std::move(name)), lambda(lambda), mu(mu), C(C) {}

  inline Real f(const Real & x) const override {
    return (this->lambda / 2. *
                (std::log(x) + std::log(this->C(0, 0) * this->C(1, 1) -
                                        Math::pow<2>(this->C(0, 1)))) +
            this->mu * (x - 1.));
  }

  inline Real f_prime(const Real & x) const override {
    AKANTU_DEBUG_ASSERT(std::abs(x) > Math::getTolerance(),
                        "x is zero (x should be the off plane right Cauchy"
                            << " measure in this function so you made a mistake"
                            << " somewhere else that lead to a zero here!!!");
    return (this->lambda / (2. * x) + this->mu);
  }

private:
  const Real & lambda;
  const Real & mu;
  const Matrix<Real> & C;
};

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void
MaterialNeohookean<dim>::computeThirdAxisDeformationOnQuad(Args && args) {
  // Neo hookean book
  auto F = Material::gradUToF<dim>(tuple::get<"grad_u"_h>(args));
  auto C = Material::rightCauchy<dim>(F);

  Math::NewtonRaphson<Real> nr(1e-5, 100);
  auto & C33 = tuple::get<"C33"_h>(args);

  C33 = nr.solve(C33_NR("Neohookean_plan_stress", this->lambda, this->mu, C),
                 C33);
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void
MaterialNeohookean<dim>::computePiolaKirchhoffOnQuad(const Matrix<Real> & E,
                                                     Matrix<Real> & S) {
  /// \f$ \sigma_{ij} = \lambda * (\nabla u)_{kk} * \delta_{ij} + \mu * (\nabla
  /// u_{ij} + \nabla u_{ji}) \f$
  S = Matrix<Real, dim, dim>::Identity() * lambda * E.trace() + 2.0 * mu * E;
}

/* -------------------------------------------------------------------------- */
template <Int dim>
inline void MaterialNeohookean<dim>::computeFirstPiolaKirchhoffOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & S, Matrix<Real> & P) {
  auto F = Material::gradUToF<dim>(grad_u);

  // first Piola-Kirchhoff is computed as the product of the deformation
  // gracient
  // tensor and the second Piola-Kirchhoff stress tensor
  P = F * S;
}

/**************************************************************************************/
/*  Computation of the potential energy for a this neo hookean material */
template <Int dim>
inline void MaterialNeohookean<dim>::computePotentialEnergyOnQuad(
    const Matrix<Real> & grad_u, Real & epot) {
  auto F = Material::gradUToF<dim>(grad_u);
  auto C = Material::rightCauchy<dim>(F);
  auto J = F.determinant();

  epot =
      0.5 * lambda * pow(log(J), 2.) + mu * (-log(J) + 0.5 * (C.trace() - dim));
}

/* -------------------------------------------------------------------------- */
template <Int dim>
template <class Args>
inline void MaterialNeohookean<dim>::computeTangentModuliOnQuad(Args && args) {
  auto & tangent = tuple::get<"tangent_moduli"_h>(args);
  // Neo hookean book
  auto cols = tangent.cols();
  auto rows = tangent.rows();
  auto F = Material::gradUToF<dim>(tuple::get<"grad_u"_h>(args));
  auto C = Material::rightCauchy<dim>(F);
  auto J = F.determinant() * sqrt(C33);

  auto && Cminus = C.inverse();

  for (Int m = 0; m < rows; m++) {
    UInt i = VoigtHelper<dim>::vec[m][0];
    UInt j = VoigtHelper<dim>::vec[m][1];
    for (Int n = 0; n < cols; n++) {
      UInt k = VoigtHelper<dim>::vec[n][0];
      UInt l = VoigtHelper<dim>::vec[n][1];

      // book belytchko
      tangent(m, n) = lambda * Cminus(i, j) * Cminus(k, l) +
                      (mu - lambda * log(J)) * (Cminus(i, k) * Cminus(j, l) +
                                                Cminus(i, l) * Cminus(k, j));
    }
  }
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
