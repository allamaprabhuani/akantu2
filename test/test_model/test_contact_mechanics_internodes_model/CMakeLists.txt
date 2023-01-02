#===============================================================================
# @file   CMakeLists.txt
#
# @author Moritz Waldleben <moritz.waldleben@gmail.com>
#
# @date creation: Thu 16 2022 
# @date last modification: Mon Jan 02 2023
#
# @brief  configuration python contact mechanics internodes test
#
#
# @section LICENSE
#
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

register_test(test_contact_detector_internodes
  SCRIPT test_contact_detector_internodes.py
  PYTHON
  FILES_TO_COPY contact2d_material.dat contact2d_plane_circle.msh
  PACKAGE python_interface contact_mechanics_internodes
  )

register_test(test_contact_mechanics_internodes_model
  SCRIPT test_contact_mechanics_internodes_model.py
  PYTHON
  FILES_TO_COPY contact2d_material.dat contact2d_plane_circle.msh contact2d_circle_circle.msh
  PACKAGE python_interface contact_mechanics_internodes
  )
