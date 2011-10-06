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
add_subdirectory("patch_tests")

#===============================================================================
set(LIST_TYPES
  triangle_3
  quadrangle_4
  tetrahedron_4
  #hexahedron_8
  )

#===============================================================================
foreach(_type ${LIST_TYPES})
  register_test(test_contact_rigid_no_friction_${_type} test_contact_rigid_no_friction.cc io_helper_tools.cc)
  set_target_properties(test_contact_rigid_no_friction_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ./data/${_type}.msh DESTINATION .)
endforeach()

#===============================================================================
file(COPY material.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)

##===============================================================================
#add_mesh(test_sliding_cube_2d_mesh sliding_cube_2d.geo 2 1)
#
#register_test(test_contact_rigid_sliding_cube_2d test_contact_rigid_sliding_cube_2d.cc)
#add_dependencies(test_contact_rigid_sliding_cube_2d test_sliding_cube_2d_mesh)
#
##===============================================================================
#add_mesh(test_sliding_cube_3d_mesh sliding_cube_3d.geo 3 1)
#
#register_test(test_contact_rigid_sliding_cube_3d test_contact_rigid_sliding_cube_3d.cc)
#add_dependencies(test_contact_rigid_sliding_cube_3d test_sliding_cube_3d_mesh)
#
##===============================================================================
#
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_2d)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_3d)

