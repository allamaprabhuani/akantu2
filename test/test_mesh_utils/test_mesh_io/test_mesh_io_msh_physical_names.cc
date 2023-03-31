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

#include "aka_common.hh"
#include "mesh.hh"
#include <iostream>
#include <sstream>

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  UInt spatialDimension(3);

  akantu::initialize(argc, argv);

  Mesh mesh(spatialDimension);

  mesh.read("./cube_physical_names.msh");
  std::stringstream sstr;

  for (auto type : mesh.elementTypes()) {
    const Array<std::string> & name_vec =
        mesh.getData<std::string>("physical_names", type);
    for (Int i(0); i < name_vec.size(); i++) {
      std::cout << "Element " << i << " (of type " << type
                << ") has physical name " << name_vec(i) << "." << std::endl;
    }
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
