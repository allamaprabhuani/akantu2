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
# @section DESCRIPTION
#
#===============================================================================


add_mesh(test_mesh plate.geo 2 1 OUTPUT plate.msh)

register_test(test_igfem_triangle_5
  SOURCES test_igfem_triangle_5.cc
  DEPENDS test_mesh
  FILES_TO_COPY material.dat
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_igfem_triangle_4
  SOURCES test_igfem_triangle_4.cc
  FILES_TO_COPY material.dat test_mesh.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )

register_test(test_interface_position
  SOURCES test_interface_position.cc
  FILES_TO_COPY material_test_interface_position.dat test_interface_position.msh
  DIRECTORIES_TO_CREATE paraview
  PACKAGE IGFEM
  )
