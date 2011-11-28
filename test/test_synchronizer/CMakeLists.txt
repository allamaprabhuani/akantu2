#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
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

if (AKANTU_SCOTCH OR AKANTU_PTSCOTCH)
  add_mesh(test_synchronizer_communication_mesh
    triangle.geo 2 1)

  register_test(test_synchronizer_communication
    test_synchronizer_communication.cc)
  add_dependencies(test_synchronizer_communication
    test_synchronizer_communication_mesh)


  register_test(test_grid_synchronizer
    test_grid_synchronizer.cc)
  add_dependencies(test_grid_synchronizer
    test_synchronizer_communication_mesh)



  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test_synchronizer_communication.sh
    ${CMAKE_CURRENT_BINARY_DIR}/test_synchronizer_communication.sh
    COPYONLY
    )

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)

  #===============================================================================
  register_test(test_dof_synchronizer test_dof_synchronizer.cc)
  copy_files(test_dof_synchronizer bar.msh)
endif()