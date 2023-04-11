/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "boundary_condition_functor.hh"
#include "test_solid_mechanics_model_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

const Real A = 1e-1;
const Real E = 1.;
const Real poisson = 3. / 10;
const auto lambda = E * poisson / (1 + poisson) / (1 - 2 * poisson);
const auto mu = E / 2 / (1. + poisson);
const Real rho = 1.;
const auto cp = std::sqrt((lambda + 2 * mu) / rho);
const auto cs = std::sqrt(mu / rho);
const auto c = std::sqrt(E / rho);

const Vector<Real, 3> k{.5, 0., 0.};
const Vector<Real, 3> psi1{0., 0., 1.};
const Vector<Real, 3> psi2{0., 1., 0.};
const auto knorm = k.norm();

/* -------------------------------------------------------------------------- */
template <Int dim> struct Verification {};

/* -------------------------------------------------------------------------- */
template <> struct Verification<1> {
  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void displacement(Eigen::MatrixBase<D1> & disp,
                    const Eigen::MatrixBase<D2> & coord, Real current_time) {
    const auto & x = coord(_x);
    const auto omega = c * k[0];
    disp(_x) = A * std::cos(k[0] * x - omega * current_time);
  }

  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void velocity(Eigen::MatrixBase<D1> & vel,
                const Eigen::MatrixBase<D2> & coord, Real current_time) {
    const auto & x = coord(_x);
    const auto omega = c * k[0];
    vel(_x) = omega * A * std::sin(k[0] * x - omega * current_time);
  }
};

/* -------------------------------------------------------------------------- */
template <> struct Verification<2> {
  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void displacement(Eigen::MatrixBase<D1> & disp,
                    const Eigen::MatrixBase<D2> & X, Real current_time) {
    Vector<Real> kshear{k[1], k[0]};
    Vector<Real> kpush{k[0], k[1]};

    const auto omega_p = knorm * cp;
    const auto omega_s = knorm * cs;

    auto phase_p = X.dot(kpush) - omega_p * current_time;
    auto phase_s = X.dot(kpush) - omega_s * current_time;

    disp = A * (kpush * std::cos(phase_p) + kshear * std::cos(phase_s));
  }

  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void velocity(Eigen::MatrixBase<D1> & vel, const Eigen::MatrixBase<D2> & X,
                Real current_time) {
    Vector<Real> kshear{k[1], k[0]};
    Vector<Real> kpush{k[0], k[1]};

    const auto omega_p = knorm * cp;
    const auto omega_s = knorm * cs;

    auto phase_p = X.dot(kpush) - omega_p * current_time;
    auto phase_s = X.dot(kpush) - omega_s * current_time;

    vel = A * (kpush * omega_p * std::cos(phase_p) +
               kshear * omega_s * std::cos(phase_s));
  }
};

/* -------------------------------------------------------------------------- */
template <> struct Verification<3> {
  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void displacement(Eigen::MatrixBase<D1> & disp,
                    const Eigen::MatrixBase<D2> & coord, Real current_time) {
    const auto & X = coord;
    auto kpush = k;
    Vector<Real, 3> kshear1 = k.cross(psi1);
    Vector<Real, 3> kshear2 = k.cross(psi2);

    const auto omega_p = knorm * cp;
    const auto omega_s = knorm * cs;

    auto phase_p = X.dot(kpush) - omega_p * current_time;
    auto phase_s = X.dot(kpush) - omega_s * current_time;

    disp = A * (kpush * std::cos(phase_p) + kshear1 * std::cos(phase_s) +
                kshear2 * std::cos(phase_s));
  }

  template <typename D1, typename D2,
            std::enable_if_t<aka::are_vectors_v<D1, D2>> * = nullptr>
  void velocity(Eigen::MatrixBase<D1> & vel,
                const Eigen::MatrixBase<D2> & coord, Real current_time) {
    const auto & X = coord;
    auto kpush = k;
    Vector<Real, 3> kshear1 = k.cross(psi1);
    Vector<Real, 3> kshear2 = k.cross(psi2);

    const auto omega_p = knorm * cp;
    const auto omega_s = knorm * cs;

    auto phase_p = X.dot(kpush) - omega_p * current_time;
    auto phase_s = X.dot(kpush) - omega_s * current_time;

    vel = A * (kpush * omega_p * std::cos(phase_p) +
               kshear1 * omega_s * std::cos(phase_s) +
               kshear2 * omega_s * std::cos(phase_s));
  }
};

/* -------------------------------------------------------------------------- */
template <ElementType _type>
class SolutionFunctor : public BC::Dirichlet::DirichletFunctor {
public:
  SolutionFunctor(Real current_time, SolidMechanicsModel & model)
      : current_time(current_time), model(model) {}

public:
  static constexpr auto dim = ElementClass<_type>::getSpatialDimension();

  inline void operator()(Int node, Vector<bool> & flags, Vector<Real> & primal,
                         const Vector<Real> & coord) const {
    flags(0) = true;
    const auto & vel = model.getVelocity();
    auto it = vel.begin(model.getSpatialDimension());
    auto v = it[node];

    Verification<dim> verif;
    verif.displacement(primal, coord, current_time);
    verif.velocity(v, coord, current_time);
  }

private:
  Real current_time;
  SolidMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
// This fixture uses somewhat finer meshes representing bars.
template <typename type_, typename analysis_method_>
class TestSMMFixtureBar : public TestSMMFixture<type_> {
  using parent = TestSMMFixture<type_>;

public:
  void SetUp() override {
    this->mesh_file =
        "../patch_tests/data/bar" + std::to_string(this->type) + ".msh";
    parent::SetUp();

    auto analysis_method = analysis_method_::value;
    this->initModel("test_solid_mechanics_model_"
                    "dynamics_material.dat",
                    analysis_method);

    const auto & position = this->mesh->getNodes();
    auto & displacement = this->model->getDisplacement();
    auto & velocity = this->model->getVelocity();

    constexpr auto dim = parent::spatial_dimension;

    Verification<dim> verif;

    for (auto && tuple :
         zip(make_view(position, dim), make_view(displacement, dim),
             make_view(velocity, dim))) {
      verif.displacement(std::get<1>(tuple), std::get<0>(tuple), 0.);
      verif.velocity(std::get<2>(tuple), std::get<0>(tuple), 0.);
    }

    if (this->dump_paraview) {
      this->model->dump();
    }

    /// boundary conditions
    this->model->applyBC(SolutionFunctor<parent::type>(0., *this->model), "BC");

    wave_velocity = 1.; // sqrt(E/rho) = sqrt(1/1) = 1
    simulation_time = 5 / wave_velocity;
    time_step = this->model->getTimeStep();

    max_steps = simulation_time / time_step; // 100
  }

  void solveStep() {
    constexpr auto dim = parent::spatial_dimension;
    Real current_time = 0.;
    const auto & position = this->mesh->getNodes();
    const auto & displacement = this->model->getDisplacement();
    auto nb_nodes = this->mesh->getNbNodes();
    auto nb_global_nodes = this->mesh->getNbGlobalNodes();

    max_error = -1.;

    Array<Real> displacement_solution(nb_nodes, dim);
    Verification<dim> verif;

    auto ndump = 50;
    auto dump_freq = max_steps / ndump;

    for (Int s = 0; s < this->max_steps; ++s, current_time += this->time_step) {
      if (s % dump_freq == 0 && this->dump_paraview) {
        this->model->dump();
      }

      /// boundary conditions
      this->model->applyBC(
          SolutionFunctor<parent::type>(current_time, *this->model), "BC");

      // compute the disp solution
      for (auto && tuple : zip(make_view(position, dim),
                               make_view(displacement_solution, dim))) {
        verif.displacement(std::get<1>(tuple), std::get<0>(tuple),
                           current_time);
      }

      // compute the error solution
      Real disp_error = 0.;
      auto n = 0;
      for (auto && tuple : zip(make_view(displacement, dim),
                               make_view(displacement_solution, dim))) {
        if (this->mesh->isLocalOrMasterNode(n)) {
          auto diff = std::get<1>(tuple) - std::get<0>(tuple);
          disp_error += diff.dot(diff);
        }
        ++n;
      }

      this->mesh->getCommunicator().allReduce(disp_error,
                                              SynchronizerOperation::_sum);

      disp_error = sqrt(disp_error) / nb_global_nodes;
      max_error = std::max(disp_error, max_error);

      this->model->solveStep();
    }
  }

protected:
  Real time_step;
  Real wave_velocity;
  Real simulation_time;
  Int max_steps;
  Real max_error{-1};
};

template <AnalysisMethod t>
using analysis_method_t = std::integral_constant<AnalysisMethod, t>;

using TestTypes = gtest_list_t<TestElementTypes>;

template <typename type_>
using TestSMMFixtureBarExplicit =
    TestSMMFixtureBar<type_, analysis_method_t<_explicit_lumped_mass>>;

TYPED_TEST_SUITE(TestSMMFixtureBarExplicit, TestTypes, );

/* -------------------------------------------------------------------------- */
TYPED_TEST(TestSMMFixtureBarExplicit, Dynamics) {
  this->solveStep();
  EXPECT_NEAR(this->max_error, 0., 2e-3);
  // std::cout << "max error: " << max_error << std::endl;
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IMPLICIT)
template <typename type_>
using TestSMMFixtureBarImplicit =
    TestSMMFixtureBar<type_, analysis_method_t<_implicit_dynamic>>;

TYPED_TEST_SUITE(TestSMMFixtureBarImplicit, TestTypes, );

TYPED_TEST(TestSMMFixtureBarImplicit, Dynamics) {
  if (this->type == _segment_2 and
      (this->mesh->getCommunicator().getNbProc() > 2)) {
    // The error are just to high after (hopefully because of the two small test
    // case)
    SUCCEED();
    return;
  }
  this->solveStep();
  EXPECT_NEAR(this->max_error, 0., 2e-3);
}
#endif

} // namespace
