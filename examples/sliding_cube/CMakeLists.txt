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
add_mesh(mesh_sliding_2d sliding_cube_2d.geo 2 1)

register_example(sliding_cube_2d sliding_cube_2d.cc)
add_dependencies(sliding_cube_2d mesh_sliding_2d )

#===============================================================================
add_mesh(mesh_sliding_3d sliding_cube_3d.geo 3 1)

register_example(sliding_cube_3d sliding_cube_3d.cc)
add_dependencies(sliding_cube_3d mesh_sliding_3d )

#===============================================================================
file(COPY material.dat DESTINATION .)
file(COPY material_elastic_caughey.dat DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_files)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_2d)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/sliding_cube_3d)
