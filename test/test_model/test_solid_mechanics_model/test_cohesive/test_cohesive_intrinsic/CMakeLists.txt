#===============================================================================
# @file   CMakeLists.txt
#
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date   Tue May 08 13:01:18 2012
#
# @brief  test for intrinsic cohesive element configuration
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

register_test(test_cohesive_intrinsic test_cohesive_intrinsic.cc)
register_test(test_cohesive_intrinsic_quadrangle test_cohesive_intrinsic_quadrangle.cc)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh.msh
  ${CMAKE_CURRENT_BINARY_DIR}/mesh.msh
  COPYONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrangle.msh
  ${CMAKE_CURRENT_BINARY_DIR}/quadrangle.msh
  COPYONLY
  )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
