#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Thu Oct 06 16:26:30 2011
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
add_mesh(mesh_hertz_2d hertz_2d.geo 2 1)

register_example(hertz_2d_no_friction hertz_2d_no_friction.cc)
add_dependencies(hertz_2d_no_friction mesh_hertz_2d )

#===============================================================================
add_mesh(mesh_hertz_3d hertz_3d.geo 3 1)

register_example(hertz_3d_no_friction hertz_3d_no_friction.cc)
add_dependencies(hertz_3d_no_friction mesh_hertz_3d )

#===============================================================================
file(COPY material.dat DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/hertz_2d_no_friction)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/hertz_3d_no_friction)
