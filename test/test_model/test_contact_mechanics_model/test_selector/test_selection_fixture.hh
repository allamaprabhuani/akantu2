/**
 * @file   test_selection_fixture.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed May 08 2019
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Generic test for selection of nodes across different surfaces
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
#include "aka_iterators.hh"
#include "communicator.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <vector>
/* -------------------------------------------------------------------------- */

using namespace akantu;

template <typename dim_>
class TestSurfaceSelecionFixture : public ::testing::Test {
public:
  static constexpr size_t spatial_dimension = dim_::value;

  void Setup() override {
    mesh = std::make_unique<Mesh>(spatial_dimension);
    mesh->read("selection_" + std::string(spatial_dimension) + "D.msh");
    selector = make_unique<PhysicalSurfaceSelector>(*mesh);
  }

  void TearDown() {
    mesh.reset(nullptr);
    selector.reset(nullptr);
  }

  bool checkNbSlaveNodes(UInt & act_slave_nodes) {
    auto & slave_list = this->selectior->getSlaveList();
    auto nb_slave_nodes = slave_list.size();
    EXCEP_EQ(nb_slave_nodes, act_slave_nodes)
  }

  bool checkNbmasterNodes(UInt & act_master_nodes) {
    auto & master_list = this->selectior->getMasterList();
    auto nb_master_nodes = master_list.size();
    EXCEP_EQ(nb_master_nodes, act_master_nodes)
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<PhysicalSurfaceSelector> selector
};

using dim_types = gtest_list_t<std::tuple<1, 2, 3>>;

TYPED_TEST_SUITE(TestSurfaceSelectionFixture, dim_types);
