/**
 * @file   test_static_memory.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 11:55:54 2010
 *
 * @brief  unit test for the StaticMemory class
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
#include "aka_static_memory.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::initialize(&argc, &argv);

  akantu::StaticMemory * st_mem = akantu::StaticMemory::getStaticMemory();

  akantu::Vector<int> & test_int = st_mem->smalloc<int>(0, "test_int", 1000, 3);

  test_int.resize(1050);

  test_int.resize(2000);

  std::cout << *st_mem << std::endl;

  st_mem->sfree(0, "test_int");

  akantu::finalize();

  exit(EXIT_SUCCESS);
}
