#===============================================================================
# @file   CMakeLists.txt
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_heat_cube3d_mesh1 cube.geo 3 1 OUTPUT cube1.msh)
add_mesh(test_heat_square_mesh1 square.geo 2 1 OUTPUT square1.msh)
#add_mesh(test_heat_bar1_mesh1 bar_structured.geo 2 1 OUTPUT bar1.msh)
add_mesh(test_heat_line_mesh line.geo 1 1 OUTPUT line.msh)
register_test(test_heat_transfer_model_cube3d
  test_heat_transfer_model_cube3d.cc)
add_dependencies(test_heat_transfer_model_cube3d test_heat_cube3d_mesh1 test_heat_square_mesh1 test_heat_line_mesh)


add_mesh(test_heat_cube_tet4 cube_tet.geo 3 1 OUTPUT cube_tet4.msh)
register_test(test_heat_transfer_model_cube3d_pbc
  test_heat_transfer_model_cube3d_pbc.cc)
add_dependencies(test_heat_transfer_model_cube3d_pbc test_heat_cube_tet4)

add_mesh(test_heat_square_tri3 square_tri.geo 2 1 OUTPUT square_tri3.msh)
register_test(test_heat_transfer_model_square2d_pbc
  test_heat_transfer_model_square2d_pbc.cc)
add_dependencies(test_heat_transfer_model_square2d_pbc test_heat_square_tri3)
register_test(test_heat_transfer_model_square2d
  test_heat_transfer_model_square2d.cc)
add_dependencies(test_heat_transfer_model_square2d test_heat_square_tri3)


register_test(test_heat_transfer_model_cube3d_istropic_conductivity
  test_heat_transfer_model_cube3d_istropic_conductivity.cc)
add_dependencies(test_heat_transfer_model_cube3d_istropic_conductivity test_heat_transfer_model_cube3d_istropic_conductivity_mesh1)

#register_test(test_heat_transfer_model_cube3d_istropic_conductivity
#  test_heat_transfer_model_cube3d_anisotropy_conductivity.cc)
#add_dependencies(test_heat_transfer_model_cube3d_anisotropy_conductivity test_heat_transfer_model_cube3d_anisotropy_conductivity_mesh1)
#
#register_test(test_heat_transfer_model_cube3d_tetra10
#  test_heat_transfer_model_cube3d_tetra10.cc)
#add_dependencies(test_heat_transfer_model_cube3d_tetra10 test_heat_cube3d_mesh2)

#===============================================================================
file(COPY material.dat DESTINATION .)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
