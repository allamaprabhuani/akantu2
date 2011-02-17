#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Wed Jan 19 12:37:24 2011
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
add_mesh(test_squares_mesh squares.geo 2 1)

register_test(test_contact_rigid_no_friction_triangle_3 test_contact_rigid_no_friction_triangle_3.cc)
add_dependencies(test_contact_rigid_no_friction_triangle_3 test_squares_mesh)

#===============================================================================
add_mesh(test_cubes_mesh cubes.geo 3 1)

register_test(test_contact_rigid_no_friction_tetrahedron_4 test_contact_rigid_no_friction_tetrahedron_4.cc)
add_dependencies(test_contact_rigid_no_friction_tetrahedron_4 test_cubes_mesh)

#===============================================================================
add_mesh(test_hertz_2d_mesh hertz_2d.geo 2 1)

register_test(test_contact_rigid_no_friction_hertz_2d test_contact_rigid_no_friction_hertz_2d.cc)
add_dependencies(test_contact_rigid_no_friction_hertz_2d test_hertz_2d_mesh)

#===============================================================================
add_mesh(test_hertz_3d_mesh hertz_3d.geo 3 1)

register_test(test_contact_rigid_no_friction_hertz_3d test_contact_rigid_no_friction_hertz_3d.cc)
add_dependencies(test_contact_rigid_no_friction_hertz_3d test_hertz_3d_mesh)

#===============================================================================
add_mesh(test_hertz_3d_full_mesh hertz_3d_full.geo 3 1)

register_test(test_contact_rigid_no_friction_hertz_3d_full test_contact_rigid_no_friction_hertz_3d_full.cc)
add_dependencies(test_contact_rigid_no_friction_hertz_3d_full test_hertz_3d_full_mesh)

#===============================================================================
add_mesh(test_force_2d_mesh force_2d.geo 2 1)

register_test(test_contact_rigid_no_friction_force_2d test_contact_rigid_no_friction_force_2d.cc)
add_dependencies(test_contact_rigid_no_friction_force_2d test_force_2d_mesh)

#===============================================================================
add_mesh(test_force_3d_mesh force_3d.geo 3 1)

register_test(test_contact_rigid_no_friction_force_3d test_contact_rigid_no_friction_force_3d.cc)
add_dependencies(test_contact_rigid_no_friction_force_3d test_force_3d_mesh)

#===============================================================================
add_mesh(test_sliding_cube_2d_mesh sliding_cube_2d.geo 2 1)

register_test(test_contact_rigid_sliding_cube_2d test_contact_rigid_sliding_cube_2d.cc)
add_dependencies(test_contact_rigid_sliding_cube_2d test_sliding_cube_2d_mesh)

#===============================================================================
add_mesh(test_sliding_cube_3d_mesh sliding_cube_3d.geo 3 1)

register_test(test_contact_rigid_sliding_cube_3d test_contact_rigid_sliding_cube_3d.cc)
add_dependencies(test_contact_rigid_sliding_cube_3d test_sliding_cube_3d_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/hertz_2d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/hertz_3d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/hertz_3d_full)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/force_2d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/force_3d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_2d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_3d)

