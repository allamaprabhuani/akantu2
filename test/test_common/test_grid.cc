/**
 * @file   test_grid.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 06 14:09:13 2012
 *
 * @brief  Test the grid object
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_grid.hh"
#include "mesh.hh"
#include "mesh_io.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  Real lower[3] = {-1.1, 0.2, -2.};
  Real upper[3] = {1.1, 1.9, 2.};

  Real spacing[3] = {0.3, 0.05, 0.5};

  RegularGrid<Real> grid(3, lower, upper, spacing);

  std::cout << grid << std::endl;

  Mesh mesh(3);

  grid.saveAsMesh(mesh);
  MeshIOMSH mesh_io;
  mesh_io.write("grid.msh", mesh);

  types::Vector<Real> position(3);

  position(0) = 0.;
  position(0) = 0.2;
  position(0) = 0.;

  grid.insert(0., position);


  akantu::finalize();

  return EXIT_SUCCESS;
}
