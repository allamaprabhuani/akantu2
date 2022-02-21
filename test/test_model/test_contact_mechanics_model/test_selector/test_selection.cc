/**
 * @file   test_selection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Generic test for surface selection
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
#include "test_selection_fixture.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

TYPED_TEST(TestSurfaceSelectionFixture, PhysicalSurface) {

  auto & actual_slave_nodes = this->mesh.getElementgroup("slave").getNodes();
  auto nb_slave_nodes = actual_slave_nodes.size();
  this->checkNbSlaveNodes(nb_slave_nodes);

  auto & actual_master_nodes = this->mesh.getElementgroup("master").getNodes();
  auto nb_master_nodes = actual_master_nodes.size();
  this->checkNbMasterNodes(nb_master_nodes);
}
