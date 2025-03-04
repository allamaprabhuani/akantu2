#===============================================================================
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
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
# @section DESCRIPTION
#
#===============================================================================


add_mesh(3d_spherical_inclusion 3d_spherical_inclusion.geo 3 2)

register_test(test_cohesive_insertion_along_physical_surfaces
  SOURCES test_cohesive_insertion_along_physical_surfaces.cc
  DEPENDS 3d_spherical_inclusion
  FILES_TO_COPY input_file.dat
  PACKAGE cohesive_element
  )

add_mesh(group_update_bar bar.geo DIM 2 ORDER 1)
register_test(test_cohesive_group_update
  SOURCES test_group_update.cc
  DEPENDS group_update_bar
  FILES_TO_COPY material.dat
  PACKAGE cohesive_element
  )
