#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.anciaux@epfl.ch>
#
# @date creation: Mon Jan 18 2016
#
# @brief  configuration for heat transfer example
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
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

add_mesh(hertz_2d hertz_2d.geo 2 1)
add_mesh(hertz_3d hertz_3d.geo 3 1)

register_example(contact_detection_2d
  SOURCES contact_detection_2d.cc
  DEPENDS hertz_2d
  FILES_TO_COPY material.dat
  )

register_example(contact_detection_3d
  SOURCES contact_detection_3d.cc
  DEPENDS hertz_3d
  FILES_TO_COPY material.dat
  )
