/**
 * @file   test_vector.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jun 29 17:32:23 2010
 *
 * @brief  test of the vector class
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
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  int def_value[3];
  def_value[0] = 10;
  def_value[1] = 20;
  def_value[2] = 30;

  akantu::Vector<int> int_vect(10, 3, def_value, "test1");

  for (unsigned int i = 5; i < int_vect.getSize(); ++i) {
    for (unsigned int j = 0; j < int_vect.getNbComponent(); ++j) {
      int_vect.values[i*int_vect.getNbComponent() + j] = def_value[j]*10;
    }
  }

  std::cerr << int_vect;

  int new_elem[3];
  new_elem[0] = 1;
  new_elem[1] = 2;
  new_elem[2] = 3;
  int_vect.push_back(new_elem);

  int_vect.push_back(200);

  int_vect.erase(0);

  std::cerr << int_vect;
  akantu::Vector<int> int_vect0(0, 3, def_value, "test2");
  std::cerr << int_vect0;
  int_vect0.push_back(new_elem);
  std::cerr << int_vect0;

  return EXIT_SUCCESS;
}
