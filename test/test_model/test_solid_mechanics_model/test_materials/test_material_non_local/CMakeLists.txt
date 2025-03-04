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


add_mesh(test_material_non_local_mesh mesh.geo 2 1 OUTPUT mesh.msh)

register_test(test_material_damage_non_local
  SOURCES test_material_damage_non_local.cc
  FILES_TO_COPY material_damage_non_local.dat mesh_section_gap.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE damage_non_local
  )


register_test(test_material_non_local
  SOURCES test_material_non_local.cc custom_non_local_test_material.cc custom_non_local_test_material.hh
  DEPENDS test_material_non_local_mesh
  FILES_TO_COPY material.dat
  PACKAGE damage_non_local
  )
