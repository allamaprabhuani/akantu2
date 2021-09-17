/**
 * @file   test_explicit_resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Thu Jan 17 2019
 *
 * @brief  Test for explicit contact resolution
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
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

const UInt spatial_dimension = 2;

int main(int argc, char * argv[]) {
  initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  // mesh.read("explicit_2d.msh");

  ContactMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  std::cout << model;

  finalize();
  return EXIT_SUCCESS;
}
