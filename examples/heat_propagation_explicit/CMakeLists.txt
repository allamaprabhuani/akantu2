#===============================================================================
# @file   CMakeLists.txt
# @author Leo
# @date   Fri Sep 29 16:46:30 2010 ...più precisi di così
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
add_mesh(double_cube_tet4 double_cube_tet.geo 3 1 OUTPUT double_cube_tet4.msh)
add_mesh(cube_tet4 cube_tet.geo 3 1 OUTPUT cube_tet4.msh)
register_example(heat_propa_pbc heat_propa_pbc.cc)
add_dependencies(heat_propa_pbc double_cube_tet4 cube_tet4)

register_example(heat_propa_parallel heat_propa_parallel.cc)
add_dependencies(heat_propa_parallel double_cube_tet4 cube_tet4)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/double_cube_tet.geo
  ${CMAKE_CURRENT_BINARY_DIR}/double_cube_tet.geo
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cube_tet.geo
  ${CMAKE_CURRENT_BINARY_DIR}/cube_tet.geo
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/heat_propa_parallel.sub
  ${CMAKE_CURRENT_BINARY_DIR}/heat_propa_parallel.sub
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/generate-meshes.sh
  ${CMAKE_CURRENT_BINARY_DIR}/generate-meshes.sh
  COPYONLY
  )

