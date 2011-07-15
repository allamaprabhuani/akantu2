# @file   CMakeLists.txt
# @author Fabian Barras <fabian.barras@epfl.ch>
# @date   Thu May 12 19:05:02 2011
#
# @brief  
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------

register_test(test_structural_mechanics_model_bernoulli_beam_2
  test_structural_mechanics_model_bernoulli_beam_2.cc)

register_test(test_structural_mechanics_model_boundary_bernoulli_beam_2
  test_structural_mechanics_model_boundary_bernoulli_beam_2.cc)

register_test(test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1
  test_structural_mechanics_model_bernoulli_beam_2_exemple_1_1.cc)


add_mesh(complicated_mesh complicated.geo 1 1)
register_test(test_structural_mechanics_model_bernoulli_beam_2_complicated
  test_structural_mechanics_model_bernoulli_beam_2_complicated.cc)
add_dependencies(test_structural_mechanics_model_bernoulli_beam_2_complicated complicated_mesh)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
