/**
 * @file   test_detection_fixture.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Dec 11 2020
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Contact detection test fixture
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "contact_mechanics_mdoel.hh"
#include "solid_mechanics_model.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TEST_DETECTION_FIXTURE_HH__
#define __AKANTU_TEST_DETECTION_FIXTURE_HH__

using namespace akantu;

template <::akantu::AnalysisMethod t>
using analysis_method_t std::integral_constant<::akantu::AnalysisMethod, t>;

class StrainIncrement : public BC::Functor {
public:
  StrainIncrement(const Matrix<Real> & strain, BC::Axis dir)
      : strain_inc(strain), dir(dir) {}

  void operator()(UInt /*node*/, Vector<bool> & flags, Vector<Real> & primal,
                  const Vector<Real> & coord) const {
    if (std::abs(coord(dir)) < 1e-8) {
      return;
    }

    flags.set(true);
    primal += strain_inc * coord;
  }

  static const BC::Functor::Type type = BC::Functor::_dirichlet;

private:
  Matrix<Real> strain_inc;
  BC::Axis dir;
};

template <typename param_> class TestCMMDFixture : public ::testing::Test {
public:
  static constexpr ElementType type_1 = std::tuple_element_t<1, param_>::value;
  static constexpr ElementType type_2 = std::tuple_element_t<2, param_>::value;

  void Setup() override {
    mesh = std::make_unique<Mesh>(this->dim);
    if (Communicator::getStaticCommunicator().whoAmI() == 0) {
      ASSERT_NO_THROW({ mesh->read(this->mesh_name); });
    }
    mesh->distribute();
  }

  void TearDown() override {
    solid.reset(nullptr);
    contact.reset(nullptr);
    mesh.reset(nullptr);
  }

  void createModel() {
    solid = std::make_unique<SolidmechanicsModel>(*mesh);
    solid->initFull(_analysis_method = this->analysis_method);

    contact = std::make_unique<ContactMechanicsModel>(*mesh);
    contact->initFull(_analysis_method = this->detection_type);
    auto && surface_selector = std::make_shared<PhysicalSurfaceSelector>(*mesh);
    contact->getContactDetector().setSurfaceSelector(surface_selector);
  }

  void setInitialCondition(const Matrix<Real> & strain) {
    for (auto && data :
         zip(make_view(this->mesh->getNodes(), this->dim),
             make_view(this->solid->getDisplacement(), this->dim))) {
      const auto & pos = std::get<0>(data);
      auto & disp = std::get<1>(data);
      disp = strain * pos;
    }
  }

  void steps(const Matrix<Real> & strain) {
    StrainIncrement functor((1. / 300) * strain, this->dim == 1 ? _x : _y);

    for (auto _ [[gnu::unused]] : arange(nb_steps)) {
      this->solid->applyBC(functor, "loading");
      this->solid->applyBC(functor, "fixed");
      this->solid->solveStep();
    }
  }

  void testImplicit() {
    this->createModel();
    auto & mat_el = this->solid->getMaterial("body");

    SCOPED_TRACE(std::to_string(this->dim) + "D - " + std::to_string(type_1) +
                 ":" + std::to_string(type_2));

    if (this->dim > 1)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "sides");
    if (this->dim > 2)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_z), "sides");

    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    Matrix<Real> strain;
    if (dim == 1) {
      strain = {{1.}};
    } else if (dim == 2) {
      strain = {{-nu, 0.}, {0., 1. - nu}};
      strain *= (1. + nu);
    } else if (dim == 3) {
      strain = {{-nu, 0., 0.}, {0., 1., 0.}, {0., 0., -nu}};
    }

    this->setInitialCondition((1 - 1e-5) * strain);
    this->steps(1e-2 * strain);

    this->contact->search();
  }

  void testExplicit() {
    this->createModel();
    auto & mat_el = this->solid->getMaterial("body");

    SCOPED_TRACE(std::to_string(this->dim) + "D - " + std::to_string(type_1) +
                 ":" + std::to_string(type_2));

    if (this->dim > 1)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "sides");
    if (this->dim > 2)
      this->model->applyBC(BC::Dirichlet::FlagOnly(_z), "sides");

    Real E = mat_el.get("E");
    Real nu = mat_el.get("nu");

    Matrix<Real> strain;
    if (dim == 1) {
      strain = {{-1.}};
    } else if (dim == 2) {
      strain = {{-nu, 0.}, {0., 1. - nu}};
      strain *= (1. + nu);
    } else if (dim == 3) {
      strain = {{-nu, 0., 0.}, {0., 1., 0.}, {0., 0., -nu}};
    }

    this->setInitialCondition((1 - 1e-5) * strain);
    this->steps(1e-2 * strain);

    this->contact->search();
  }

  bool checkGap() {}

  bool checkNormal() {}

  bool checkCovariantBasis() {}

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> solid;
  std::unique_ptr<ContactMechanicsModel> contact;

  std::string mesh_name{
      std::to_string(detection_type) + std::to_string(type_1) +
      (type_1 == type_2 ? "" : std::to_string(type_2)) + ".msh"};

  AnalysisMethod analysis_method;
  DetectionType detection_type;
};

/* -------------------------------------------------------------------------- */

using element_types = gtest_list_t<std::tuple<
    std::tuple<_element_type_segment_2, _element_type_segement_2>,
    std::tuple<_element_type_triangle_3, _element_type_triangle_3>,
    std::tuple<_element_type_triangle_3, _element_type_quadrangle_4>,
    std::tuple<_element_type_quadrangle_4, _element_type_quadrange_4>,
    std::tuple<_element_type_tetrahedron_6, _element_type_tetrahedron_6>>>;

TYPED_TEST_SUITE(TestCMMDFixture, detection_types)

#endif
