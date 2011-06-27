/**
 * @file   test_mesh_io_msh.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 14 14:27:11 2010
 *
 * @brief  unit test for the MeshIOMSH class
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
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */


int main(int argc, char *argv[]) {
  akantu::initialize(&argc,&argv);
  akantu::MeshIOMSH mesh_io;
  akantu::Mesh mesh(3);

  mesh_io.read("./cube.msh", mesh);

  std::cout << mesh << std::endl;

  akantu::finalize();
  return EXIT_SUCCESS;
}
