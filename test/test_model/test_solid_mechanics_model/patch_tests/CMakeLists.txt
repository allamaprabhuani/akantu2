#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Wed Feb 16 09:46:59 2010
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
add_mesh(test_weight_hexa_mesh colone_hexa.geo 3 1)
add_mesh(test_weight_tetra_mesh colone_tetra.geo 3 1)

register_test(test_weight colone_weight.cc)
add_dependencies(test_weight test_weight_hexa_mesh test_weight_tetra_mesh)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material_colone.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material_colone.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/test_weight_hexa)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/test_weight_tetra)
#===============================================================================

set(LIST_TYPES
  segment_2
  segment_3
  triangle_3
  triangle_6
  quadrangle_4
  quadrangle_8
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  )

add_subdirectory(explicit)

if(AKANTU_MUMPS)
  add_subdirectory(implicit)
endif()
