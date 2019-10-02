#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Tue Dec 18 2018
# @date last modification: Tue Dec 18 2018
#
# @brief  configuration for ContactMechanicsModel tests
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

add_akantu_test(test_detector "detection_test")
add_akantu_test(test_resolution "resolution_test")
add_akantu_test(test_coupler "coupler_test")

add_mesh (hertz_2d hertz_2d.geo 2 1)

register_test(test_hertz_2d
  SOURCES test_hertz_2d.cc
  DEPENDS hertz_2d
  FILES_TO_COPY material.dat
  PACKAGE contact_mechanics
  )

