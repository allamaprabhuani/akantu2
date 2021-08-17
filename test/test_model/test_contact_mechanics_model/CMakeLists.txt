#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Sun Jun 06 2021
#
# @brief  configuration for ContactMechanicsModel tests
#
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
# @section DESCRIPTION
#
#===============================================================================


add_akantu_test(test_selector "selector_test")

add_mesh(flat_on_flat flat_on_flat.geo 2 1)
register_test(test_explicit_dynamic
  SOURCES test_explicit_dynamic.cc
  DEPENDS flat_on_flat
  FILES_TO_COPY material.dat
  PACKAGE contact_mechanics
  )


add_mesh(sliding-block-2D sliding-block-2D.geo 2 1)
register_test(test_explicit_friction
  SOURCES test_explicit_friction.cc
  DEPENDS sliding-block-2D
  FILES_TO_COPY material-friction.dat
  PACKAGE contact_mechanics
  )
