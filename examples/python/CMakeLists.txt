#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Jan 22 2016
# @date last modification: Tue Mar 30 2021
#
# @brief  Main CMakeLists for python examples
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


add_subdirectory(custom-material)
add_subdirectory(dynamics)
add_subdirectory(eigen_modes)
add_subdirectory(plate-hole)
add_subdirectory(stiffness_matrix)

add_example(structural_mechanics "structural mechanics example in python"
  PACKAGE structural_mechanics)
add_example(fragmentation "example of fragmentation in python"
  PACKAGE cohesive_element)
add_example(phase-field "phase-field example in python"
  PACKAGE phase_field)
add_example(cohesive "cohesive element examples in python"
  PACKAGE cohesive_element)
add_example(contact-mechanics "contact mechanics example in python"
  PACKAGE contact_mechanics)
package_add_files_to_package(
  examples/python/README.rst
  )
