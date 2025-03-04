#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_mesh(single_agg energy_study.geo 2 1)

register_test(test_material_damage_iterative
  SOURCES test_material_damage_iterative.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials implicit
  )

register_test(test_material_damage_iterative_parallel
  SOURCES test_material_damage_iterative.cc
  FILES_TO_COPY material.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials implicit
  )

register_test(test_material_damage_iterative_non_local_serial
  SOURCES test_material_damage_iterative_non_local_serial.cc
  FILES_TO_COPY material_non_local.dat plate.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials_non_local implicit
  )

register_test(test_material_damage_iterative_non_local_parallel
  SOURCES test_material_damage_iterative_non_local_parallel.cc
  FILES_TO_COPY two_materials.dat one_circular_inclusion.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials_non_local implicit
  )

register_test(test_material_iterative_stiffness_reduction
  SOURCES test_material_iterative_stiffness_reduction.cc
  FILES_TO_COPY material_stiffness_reduction.dat two_elements.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

# register_test(test_material_iterative_stiffness_reduction_energy
#   SOURCES test_material_iterative_stiffness_reduction_energy.cc
#   DEPENDS single_agg
#   FILES_TO_COPY concrete.dat
#   DIRECTORIES_TO_CREATE paraview
#   PACKAGE extra_materials
#   )


