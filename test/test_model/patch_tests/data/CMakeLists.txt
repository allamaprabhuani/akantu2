#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


add_mesh(bar_segment_2 bar_segment.geo
  DIM 1 ORDER 1 OUTPUT bar_segment_2.msh)
add_mesh(bar_segment_3 bar_segment.geo
  DIM 1 ORDER 2 OUTPUT bar_segment_3.msh)
add_mesh(bar_triangle_3 bar_triangle.geo
  DIM 2 ORDER 1 OUTPUT bar_triangle_3.msh)
add_mesh(bar_triangle_6 bar_triangle.geo
  DIM 2 ORDER 2 OUTPUT bar_triangle_6.msh)
add_mesh(bar_quadrangle_4 bar_quadrangle.geo
  DIM 2 ORDER 1 OUTPUT bar_quadrangle_4.msh)
add_mesh(bar_quadrangle_8 bar_quadrangle.geo
  DIM 2 ORDER 2 OUTPUT bar_quadrangle_8.msh)
add_mesh(bar_tetrahedron_4 bar_tetrahedron.geo
  DIM 3 ORDER 1 OUTPUT  bar_tetrahedron_4.msh)
add_mesh(bar_tetrahedron_10 bar_tetrahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_tetrahedron_10.msh)
add_mesh(bar_hexahedron_8  bar_hexahedron.geo
  DIM 3 ORDER 1 OUTPUT bar_hexahedron_8.msh)
add_mesh(bar_hexahedron_20 bar_hexahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_hexahedron_20.msh)
add_mesh(bar_pentahedron_6 bar_pentahedron.geo
  DIM 3 ORDER 1 OUTPUT bar_pentahedron_6.msh)
add_mesh(bar_pentahedron_15 bar_pentahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_pentahedron_15.msh)

#===============================================================================
# Meshes
#===============================================================================
set(_patch_tests_meshes)
set(_patch_tests_bar_meshes)
package_get_variable(ET_ELEMENT_TYPES core _element_list)

set(_element_list
  _hexahedron_8
  _hexahedron_20
  _pentahedron_6
  _pentahedron_15
  _quadrangle_4
  _quadrangle_8
  _segment_2
  _segment_3
  _tetrahedron_4
  _tetrahedron_10
  _triangle_3
  _triangle_6
  )

foreach(_et ${_element_list})
  list(APPEND _patch_tests_meshes ${CMAKE_CURRENT_SOURCE_DIR}/${_et}.msh)
  list(APPEND _patch_tests_bar_meshes bar${_et})
endforeach()

set(PATCH_TEST_MESHES ${_patch_tests_meshes}
  CACHE INTERNAL "List of mesh files")

set(PATCH_TEST_BAR_MESHES ${_patch_tests_bar_meshes}
  CACHE INTERNAL "List of bar mesh files")
