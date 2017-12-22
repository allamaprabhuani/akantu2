/**
 * @file   test_fe_engine_precomputation.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Jul 13 2015
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "pybind11_akantu.hh"
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
/* -------------------------------------------------------------------------- */
using namespace akantu;

namespace py = pybind11;
using namespace py::literals;

template <typename type_>
class TestFEMPyFixture : public TestFEMFixture<type_> {
  using parent = TestFEMFixture<type_>;

public:
  void SetUp() override {
    parent::SetUp();

    const auto & connectivities = this->mesh->getConnectivity(this->type);
    const auto & nodes = this->mesh->getNodes().begin(this->dim);
    coordinates = std::make_unique<Array<Real>>(
        connectivities.size(), connectivities.getNbComponent() * this->dim);

    for (auto && tuple :
         zip(make_view(connectivities, connectivities.getNbComponent()),
             make_view(*coordinates, this->dim,
                       connectivities.getNbComponent()))) {
      const auto & conn = std::get<0>(tuple);
      const auto & X = std::get<1>(tuple);
      for (auto s : arange(conn.size())) {
        Vector<Real>(X(s)) = Vector<Real>(nodes[conn(s)]);
      }
    }
  }

  void TearDown() override {
    parent::TearDown();
    coordinates.reset(nullptr);
  }

protected:
  std::unique_ptr<Array<Real>> coordinates;
};

TYPED_TEST_CASE(TestFEMPyFixture, types);

TYPED_TEST(TestFEMPyFixture, DISABLED_Precompute) {
  const auto & N = this->fem->getShapeFunctions().getShapes(this->type);
  const auto & B =
      this->fem->getShapeFunctions().getShapesDerivatives(this->type);
  const auto & j = this->fem->getIntegrator().getJacobians(this->type);

  Array<Real> ref_N(this->nb_quadrature_points_total, N.getNbComponent());
  Array<Real> ref_B(this->nb_quadrature_points_total, B.getNbComponent());
  Array<Real> ref_j(this->nb_quadrature_points_total, j.getNbComponent());

  py::module py_engine = py::module::import("py_engine");
  auto py_shape = py_engine.attr("Shapes")(py::str(aka::to_string(this->type)));
  auto kwargs = py::dict("N"_a = make_proxy(ref_N), "B"_a = make_proxy(ref_B),
                         "j"_a = make_proxy(ref_j),
                         "X"_a = make_proxy(*this->coordinates));

  auto ret = py_shape.attr("precompute")(**kwargs);

  auto check = [&](auto & ref_A, auto & A) {
    for (auto && n :
         zip(make_view(ref_A,
                       this->nb_quadrature_points_total / this->nb_element,
                       ref_A.getNbComponent()),
             make_view(A, this->nb_quadrature_points_total / this->nb_element,
                       A.getNbComponent()))) {
      auto diff = (std::get<0>(n) - std::get<1>(n)).template norm<L_inf>();
      //std::cout << "||" << std::get<0>(n) << " - " << std::get<1>(n) << "|| = " << diff << std::endl;
      EXPECT_NEAR(0., diff, 1e-14);
    }
  };
  check(ref_N, N);
  check(ref_B, B);
  check(ref_j, j);
}
