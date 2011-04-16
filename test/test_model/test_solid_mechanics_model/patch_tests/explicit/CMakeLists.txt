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

register_test(test_check_stress test_check_stress.cc)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material_check_stress.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material_check_stress.dat
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/square_check_stress.msh
  ${CMAKE_CURRENT_BINARY_DIR}/square_check_stress.msh
  COPYONLY
  )

#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/check_stress_test)
