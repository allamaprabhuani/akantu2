#===============================================================================
# @file   CMakeLists.txt
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
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

add_mesh(test_mesh_partitionate_mesh
  triangle.geo 2 2)

register_test(test_mesh_partitionate_scotch
  test_mesh_partitionate_scotch.cc)
add_dependencies(test_mesh_partitionate_scotch
  test_mesh_partitionate_mesh)

register_test(test_mesh_partitionate_scotch_advanced
  test_mesh_partitionate_scotch_advanced.cc)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/squares_L.msh
  ${CMAKE_CURRENT_BINARY_DIR}/squares_L.msh
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/squares_H.msh
  ${CMAKE_CURRENT_BINARY_DIR}/squares_H.msh
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
