#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Sun Sep 12 23:37:43 2010
#
# @brief  configuration for synchronizer tests
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

add_mesh(test_synchronizer_communication_mesh
  triangle.geo 2 1)

register_test(test_synchronizer_communication
  SOURCES test_synchronizer_communication.cc
  DEPENDENCIES test_synchronizer_communication_mesh
  )

add_executable(test_grid_synchronizer_check_neighbors test_grid_synchronizer_check_neighbors.cc)
target_link_libraries(test_grid_synchronizer_check_neighbors akantu)

register_test(test_grid_synchronizer
  SOURCES test_grid_synchronizer.cc
  DEPENDENCIES test_synchronizer_communication_mesh test_grid_synchronizer_check_neighbors
  )

register_test(test_dof_synchronizer
  SOURCES test_dof_synchronizer.cc
  FILES_TO_COPY bar.msh)
