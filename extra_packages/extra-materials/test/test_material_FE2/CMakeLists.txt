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
##===============================================================================


register_test(test_material_FE2
  SOURCES test_material_FE2.cc
  FILES_TO_COPY material.dat periodic_plate.msh square.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

register_test(test_eigenstrain_homogenization
  SOURCES test_eigenstrain_homogenization.cc
  FILES_TO_COPY mesoscale_materials.dat one_inclusion.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE extra_materials
  )

add_mesh(homogenized_plate homogenized_plate.geo 2 1)
register_test(test_elastic_homogenization
  SOURCES test_elastic_homogenization.cc
  DEPENDS homogenized_plate
  FILES_TO_COPY material_orthotropic.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit extra_materials
  )

add_mesh(periodic_plate periodic_plate.geo 2 1)
register_test(test_node_selection
  SOURCES  test_node_selection.cc 
  DEPENDS periodic_plate
  FILES_TO_COPY material_test_boundary.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE implicit extra_materials
  )
