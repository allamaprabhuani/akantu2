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


add_mesh(test_msh_cube cube.geo 3 1)
add_mesh(test_msh_cube_physical_names cube_physical_names.geo 3 1)

register_test(test_mesh_io_msh
  SOURCES test_mesh_io_msh.cc
  DEPENDS test_msh_cube
  DIRECTORIES_TO_CREATE paraview
  PACKAGE core
  )

register_test(test_mesh_io_msh_physical_names
  SOURCES test_mesh_io_msh_physical_names.cc
  DEPENDS test_msh_cube_physical_names
  PACKAGE core
  )

#===============================================================================
#register_test(test_mesh_io_diana test_mesh_io_diana.cc)
#copy_files(test_mesh_io_diana dam.dat)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview/mesh_io_diana)

