/**
 * @file   test_detection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Generic test for detection between different element types
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
#include "aka_iterators.hh"
#include "communicator.hh"
#include "test_detection_fixture.hh"
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestCMMDFixture, Implicit) {

  this->analysis_method = _static;
  this->detection_type = _implicit_contact;

  this->testImplicit();

  this->checkGap();
}

TYPED_TEST(TestCMMDFixture, Explicit) {

  this->analysis_method = _static;
  this->detection_type = _explicit_contact;

  this->testExplicit();

  this->checkGap();
}
