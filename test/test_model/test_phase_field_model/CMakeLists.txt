#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Jun 08 2021
#
# @brief  test for the phase field model
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
#===============================================================================


register_test(test_phasefield_selector
  SOURCES test_phasefield_selector.cc
  FILES_TO_COPY phasefield_selector.dat phasefield_selector.msh
  PACKAGE phase_field
  )

register_test(test_phase_solid_coupling
  SOURCES test_phase_solid_coupling.cc
  FILES_TO_COPY material_coupling.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_phase_solid_explicit
  SOURCES test_phase_solid_explicit.cc
  FILES_TO_COPY material_coupling.dat test_one_element.msh
  PACKAGE phase_field
  )

register_test(test_multi_material
  SOURCES test_multi_material.cc
  FILES_TO_COPY material_multiple.dat test_two_element.msh
  PACKAGE phase_field
  )
