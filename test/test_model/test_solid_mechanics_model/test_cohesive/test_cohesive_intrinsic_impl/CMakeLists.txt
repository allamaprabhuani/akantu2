#===============================================================================
# @file   CMakeLists.txt
#
# @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
# @author Marco Vocialta <marco.vocialta@epfl.ch>
#
# @date   Mon Jul 09 14:13:56 2012
#
# @brief  configuration for intrinsic implicit cohesive elements
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

register_test(test_cohesive_intrinsic_impl test_cohesive_intrinsic_impl.cc)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )


file(COPY material.dat  DESTINATION .)
file(COPY simple.msh DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(COPY simple.msh DESTINATION .)