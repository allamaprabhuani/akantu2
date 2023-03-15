/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

const Int spatial_dimension = 2;

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
