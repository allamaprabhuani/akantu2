#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Wed Dec 18 2018
# @date last modification: Wed Apr 29 2019
#
# @brief  configuration for contact detection tests
#
# @section LICENSE
#
# Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh (node node.geo 2 1)
add_mesh (hertz_explicit hertz_explicit.geo 2 1)

register_test(test_node_selector
  SOURCES test_node_selector.cc
  DEPENDS node
  FILES_TO_COPY options.dat
  PACKAGE contact_mechanics
  )

register_test(test_explicit_detection
  SOURCES test_explicit_detection.cc
  DEPENDS hertz_explicit
  FILES_TO_COPY options_explicit.dat
  PACKAGE contact_mechanics
  )

