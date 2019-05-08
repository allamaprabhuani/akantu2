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

add_mesh (implicit_2d implicit_2d.geo 2 1)
add_mesh (spheres_2d  spheres_2d.geo  2 1)
add_mesh (implicit_3d implicit_3d.geo 3 1)


register_test(test_implicit_3d
  SOURCES test_detection_implicit_3d.cc
  DEPENDS implicit_3d
  FILES_TO_COPY options.dat
  PACKAGE contact_mechanics
  )

register_test(test_implicit_2d
  SOURCES test_detection_implicit_2d.cc
  DEPENDS implicit_2d
  FILES_TO_COPY options.dat
  PACKAGE contact_mechanics
  )

register_test(test_spheres_2d
  SOURCES test_detection_implicit_2d.cc
  DEPENDS spheres_2d
  FILES_TO_COPY options.dat
  PACKAGE contact_mechanics
  )

