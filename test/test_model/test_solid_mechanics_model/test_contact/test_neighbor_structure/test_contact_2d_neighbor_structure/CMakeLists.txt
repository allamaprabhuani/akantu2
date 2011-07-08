#===============================================================================
# @file   CMakeLists.txt
# @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
# @date   Thu Dec  9 10:19:20 2010
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
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_contact_2d_neighbor_structure_mesh squares.geo 2 1)

register_test(test_contact_2d_neighbor_structure test_contact_2d_neighbor_structure.cc)
add_dependencies(test_contact_2d_neighbor_structure test_contact_2d_neighbor_structure_mesh)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(COPY materials.dat DESTINATION .)
