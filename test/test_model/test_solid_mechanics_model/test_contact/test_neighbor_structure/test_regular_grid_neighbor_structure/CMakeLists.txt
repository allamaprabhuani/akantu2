#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Tue Oct 26 16:40:24 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
add_mesh(test_regular_grid_2d_mesh squares.geo 2 1)

register_test(test_regular_grid_triangle_3 test_regular_grid_triangle_3.cc)
add_dependencies(test_regular_grid_triangle_3 test_regular_grid_2d_mesh)

#===============================================================================
add_mesh(test_regular_grid_3d_mesh cubes.geo 3 1)

register_test(test_regular_grid_tetrahedron_4 test_regular_grid_tetrahedron_4.cc)
add_dependencies(test_regular_grid_tetrahedron_4 test_regular_grid_3d_mesh)

#===============================================================================
#add_mesh(test_regular_grid_3d_mesh cubes.geo 3 1)

register_test(test_regular_grid_tetrahedron_4_nodes test_regular_grid_tetrahedron_4_nodes.cc)
add_dependencies(test_regular_grid_tetrahedron_4_nodes test_regular_grid_3d_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
