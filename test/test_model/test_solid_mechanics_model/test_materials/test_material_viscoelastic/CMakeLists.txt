#===============================================================================
# Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================


add_mesh(test_material_standard_linear_solid_deviatoric_relaxation_mesh
  test_material_standard_linear_solid_deviatoric_relaxation.geo 2 1)

register_test(test_material_standard_linear_solid_deviatoric_relaxation
  SOURCES test_material_standard_linear_solid_deviatoric_relaxation.cc
  DEPENDS test_material_standard_linear_solid_deviatoric_relaxation_mesh
  FILES_TO_COPY material_standard_linear_solid_deviatoric_relaxation.dat
  PACKAGE core
  )

register_test(test_material_standard_linear_solid_deviatoric_relaxation_tension
  SOURCES test_material_standard_linear_solid_deviatoric_relaxation_tension.cc
  DEPENDS test_material_standard_linear_solid_deviatoric_relaxation_mesh
  FILES_TO_COPY material_standard_linear_solid_deviatoric_relaxation.dat
  PACKAGE core
  )
