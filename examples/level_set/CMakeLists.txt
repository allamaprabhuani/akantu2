#===============================================================================
# @file   CMakeLists.txt
#
# @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
#
# @date   Fri Feb 24 03:42:39 2012
#
# @brief  configuration level set test
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
register_example(level_set level_set.cc)

add_mesh(level_set_mesh square_level_set.geo 2 2)
add_dependencies(level_set level_set_mesh)

#===============================================================================
#register_example(implicit_dynamic implicit_dynamic.cc)

#add_mesh(implicit_dynamic_mesh beam.geo 2 2)
#add_dependencies(implicit_dynamic implicit_dynamic_mesh)

#===============================================================================
#configure_file(
#  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
#  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
#  COPYONLY
#  )

#configure_file(
#  ${CMAKE_CURRENT_SOURCE_DIR}/material_dynamic.dat
#  ${CMAKE_CURRENT_BINARY_DIR}/material_dynamic.dat
#  COPYONLY
#  )


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
