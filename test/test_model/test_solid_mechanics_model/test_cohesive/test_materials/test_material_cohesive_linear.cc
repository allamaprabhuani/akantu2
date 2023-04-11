/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_material_cohesive_fixture.hh"
/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear.hh"
/* -------------------------------------------------------------------------- */

template <Int dim>
struct TestMaterialCohesiveLinear
    : public TestMaterialCohesive<MaterialCohesiveLinear, dim> {

  TestMaterialCohesiveLinear(SolidMechanicsModel & model)
      : TestMaterialCohesive<MaterialCohesiveLinear, dim>(model) {}

  void SetUp() override {
    this->is_extrinsic = true;

    this->beta = 2.;
    this->kappa = 2;

    this->G_c = 10.;
    this->sigma_c_ = 1e6;
    this->penalty = 1e11;

    this->delta_c_ = 2. * this->G_c / this->sigma_c_;
  }

  void resetInternal() override {
    // normal_opening_ = Vector<Real>::Zero(dim);
    // tangential_opening_ = Vector<Real>::Zero(dim);
    contact_traction_ = Vector<Real>::Zero(dim);
    contact_opening_ = Vector<Real>::Zero(dim);
  }

  void computeTractions(Array<Real> & openings, const Vector<Real> & normal,
                        Array<Real> & tractions) override {
    for (auto && [opening, traction] :
         zip(make_view(openings, dim), make_view(tractions, dim))) {

      auto && args = tuple::make_named_tuple(
          "normal"_n = normal, "opening"_n = opening, "traction"_n = traction,
          "contact_opening"_n = contact_opening_,
          "contact_traction"_n = contact_traction_, "delta_max"_n = delta_max_,
          "sigma_c"_n = this->sigma_c_, "damage"_n = damage_,
          "delta_c"_n = this->delta_c_,
          "insertion_stress"_n = this->insertion_stress_);

      this->computeTractionOnQuad(args);

      opening += contact_opening_;
      traction += contact_traction_;
    }
  }

  Real delta(const Vector<Real> & opening, const Vector<Real> & normal) {
    auto beta = this->beta;
    auto kappa = this->kappa;

    auto normal_opening = opening.dot(normal) * normal;
    auto tangential_opening = opening - normal_opening;

    return std::sqrt(std::pow(normal_opening.norm(), 2) +
                     std::pow(tangential_opening.norm() * beta / kappa, 2));
  }

  Vector<Real> traction(const Vector<Real> & opening,
                        const Vector<Real> & normal) {
    auto delta_c = this->delta_c_;
    auto sigma_c = this->sigma_c_;
    auto beta = this->beta;
    auto kappa = this->kappa;

    auto normal_opening = opening.dot(normal) * normal;
    auto tangential_opening = opening - normal_opening;

    auto delta_ = this->delta(opening, normal);
    if (delta_ < 1e-16) {
      return this->insertion_stress_;
    }

    if (opening.dot(normal) / delta_c < -Math::getTolerance()) {
      ADD_FAILURE() << "This is contact";
      return Vector<Real>::Zero(dim);
    }

    auto T = sigma_c * (delta_c - delta_) / delta_c / delta_ *
             (normal_opening + tangential_opening * beta * beta / kappa);

    return T;
  }

  Vector<Real> tractionModeI(const Vector<Real> & opening,
                             const Vector<Real> & normal) {
    return traction(opening, normal);
  }

  Vector<Real> tractionModeII(const Vector<Real> & opening,
                              const Vector<Real> & normal) {
    return traction(opening, normal);
  }

public:
  Real delta_c_{0.};

  Real delta_max_{0.};
  // Real normal_opening_norm_{0.};
  // Real tangential_opening_norm{0.};
  Real damage_{0.};
  bool penetration_{false};

  Real etot{0.};
  Real edis{0.};

  Vector<Real> normal_opening_;
  Vector<Real> tangential_opening_;

  Vector<Real> contact_traction_;
  Vector<Real> contact_opening_;
};

template <typename dim_>
using TestMaterialCohesiveLinearFixture =
    TestMaterialCohesiveFixture<TestMaterialCohesiveLinear, dim_>;

using coh_types = gtest_list_t<TestAllDimensions>;

TYPED_TEST_SUITE(TestMaterialCohesiveLinearFixture, coh_types, );

TYPED_TEST(TestMaterialCohesiveLinearFixture, ModeI) {
  this->checkModeI(this->material->delta_c_, this->material->get("G_c"));

  Real G_c = this->material->get("G_c");
  EXPECT_NEAR(G_c, this->dissipated(), 1e-6);
}

TYPED_TEST(TestMaterialCohesiveLinearFixture, ModeII) {
  this->checkModeII(this->material->delta_c_);

  if (this->dim != 1) {
    Real G_c = this->material->get("G_c");
    Real beta = this->material->get("beta");
    Real dis = beta * G_c;
    EXPECT_NEAR(dis, this->dissipated(), 1e-6);
  }
}

TYPED_TEST(TestMaterialCohesiveLinearFixture, Cycles) {
  auto delta_c = this->material->delta_c_;
  auto sigma_c = this->material->sigma_c_;

  this->material->insertion_stress_ = this->normal * sigma_c;

  this->addOpening(this->normal, 0, 0.1 * delta_c, 100);
  this->addOpening(this->normal, 0.1 * delta_c, 0., 100);
  this->addOpening(this->normal, 0., 0.5 * delta_c, 100);
  this->addOpening(this->normal, 0.5 * delta_c, -1.e-5, 100);
  this->addOpening(this->normal, -1.e-5, 0.9 * delta_c, 100);
  this->addOpening(this->normal, 0.9 * delta_c, 0., 100);
  this->addOpening(this->normal, 0., delta_c, 100);

  this->material->computeTractions(*this->openings, this->normal,
                                   *this->tractions);

  Real G_c = this->material->get("G_c");
  EXPECT_NEAR(G_c, this->dissipated(), 2e-3); // due to contact dissipation at 0
  this->output_csv();
}
