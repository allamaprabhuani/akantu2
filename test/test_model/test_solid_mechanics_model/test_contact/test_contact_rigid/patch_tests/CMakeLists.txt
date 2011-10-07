#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Wed May 25 12:37:24 2011
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
add_mesh(test_single_spring_mesh single_triangle.geo 2 1)

register_test(test_single_spring_friction test_single_spring_friction.cc)
add_dependencies(test_single_spring_friction test_single_spring_mesh)

register_test(test_single_spring_friction_w_visco test_single_spring_friction_w_visco.cc)
add_dependencies(test_single_spring_friction_w_visco test_single_spring_mesh)

#register_test(test_single_spring_friction_vel_weak_exp test_single_spring_friction_vel_weak_exp.cc)
#add_dependencies(test_single_spring_friction_vel_weak_exp test_single_spring_mesh)

#===============================================================================
file(COPY material.dat DESTINATION .)
file(COPY material_elastic_caughey.dat DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files/single_spring_friction)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files/single_spring_friction_w_visco)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files/single_spring_friction_vel_weak_exp)
if(AKANTU_IOHELPER_ON)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/single_spring_friction)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/single_spring_friction_w_visco)
  #file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/single_spring_friction_vel_weak_exp)
endif()
