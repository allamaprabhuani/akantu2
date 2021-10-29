/**
 * @file   test_cohesive.cc
 *
 * @author Tobias Brink <tobias.brink@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jan 10 2018
 * @date last modification:  Sun Dec 30 2018
 *
 * @brief  Generic test for cohesive elements
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_cohesive_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestSMMCFixture, ExtrinsicModeI) {
  if (this->mesh->getCommunicator().getNbProc() > 1 and this->dim == 1) {
    SUCCEED();
    return;
  }
  getStaticParser().parse("material_0.dat");
  this->is_extrinsic = true;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeI();
  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  //  if (this->dim != 3)
  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, ExtrinsicModeIFiniteDef) {
  if (this->dim == 1) {
    SUCCEED();
    return;
  }
  getStaticParser().parse("material_0_finite_def.dat");
  this->is_extrinsic = true;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeI();
  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  //  if (this->dim != 3)
  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, ExtrinsicModeII) {
  if (this->mesh->getCommunicator().getNbProc() > 1 and this->dim == 1) {
    SUCCEED();
    return;
  }
  getStaticParser().parse("material_0.dat");
  this->is_extrinsic = true;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeII();
  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  // if (this->dim != 3)
  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, IntrinsicModeI) {
  if (this->mesh->getCommunicator().getNbProc() > 1 and this->dim == 1) {
    SUCCEED();
    return;
  }
  getStaticParser().parse("material_1.dat");
  this->is_extrinsic = false;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeI();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  // if (this->dim != 3)
  this->checkDissipated(G_c);
}

TYPED_TEST(TestSMMCFixture, IntrinsicModeII) {
  if (this->mesh->getCommunicator().getNbProc() > 1 and this->dim == 1) {
    SUCCEED();
    return;
  }

  getStaticParser().parse("material_1.dat");
  this->is_extrinsic = false;
  this->analysis_method = _explicit_lumped_mass;

  this->testModeII();

  this->checkInsertion();

  auto & mat_co = this->model->getMaterial("insertion");
  Real G_c = mat_co.get("G_c");

  // if (this->dim != 3)
  this->checkDissipated(G_c);
}
