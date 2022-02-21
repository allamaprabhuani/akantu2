/**
 * @file   test_resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed May 08 2019
 * @date last modification: Wed May 08 2019
 *
 * @brief  Test the resolution class common function
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
#include "resolution.hh"
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.hh>
#include <random>
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */

class TestResolutionFixture : public ::testing::Test {
public:
  void SetUp() override {
    mesh = std::make_unique<Mesh>(spatial_dimension);
    model = std::make_unique<Model>(*mesh);
    resolution = std::make_unique<Resolution>(*model, "resolution");
  }

  void TearDown() override {
    resolution.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<ContactMechanicsModel> model;
  std::unique_ptr<Resolution> resolution;
};

TYPED_TEST(TestResolutionFixture, TestComputeN) {

  Vector<Real> shapes(nb_nodes_master);
  Vector<Real> projection(spatial_diemnsion - 1);

#define GET_SHAPES_NATURAL(type)                                               \
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

  Vector<Real> n(conn.size() * spatial_dimension);
  resolution->computeN(n, shapes, normal);
}

TYPED_TEST(TestResolutionFixture, TestComputeDalpha) {}

TYPED_TEST(TestResolutionFixture, TestComputeNalpha) {}
