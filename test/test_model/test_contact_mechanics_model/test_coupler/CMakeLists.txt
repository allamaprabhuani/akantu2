#===============================================================================
# @file   CMakeLists.txt
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri May 24 2019
#
# @brief  configuration for testing coupler class for contact mechanics model
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


add_mesh (coupling coupling.geo 2 1)
add_mesh (flat_on_flat flat_on_flat.geo 2 1)


register_test( test_contact_coupling
  SOURCES test_contact_coupling.cc
  DEPENDS coupling
  FILES_TO_COPY material.dat
  PACKAGE contact_mechanics model_couplers
  )


register_test( test_coupled_stiffness
  SOURCES test_coupled_stiffness.cc
  DEPENDS flat_on_flat
  FILES_TO_COPY material_stiffness.dat
  PACKAGE contact_mechanics model_couplers
  )
