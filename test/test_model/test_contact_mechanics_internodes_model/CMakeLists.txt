#===============================================================================
# @file   CMakeLists.txt
#
# @author Moritz Waldleben <moritz.waldleben@gmail.com>
#
# @date creation: Thu 16 2022 
# @date last modification: Thue 16 2022 
#
# @brief  configuration python contact mechanics internodes test
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

add_mesh(contact_mesh prototype_internodes/contact.geo DIM 2 ORDER 1 OUTPUT contact.msh)

register_test(test_internodes
  SCRIPT test_internodes.py
  PYTHON
  FILES_TO_COPY prototype_internodes/material.dat
  DEPENDS contact_mesh
  PACKAGE python_interface contact_mechanics_internodes_model
  )
