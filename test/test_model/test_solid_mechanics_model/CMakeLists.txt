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

add_akantu_test(test_materials "test_materials")
add_akantu_test(patch_tests "patch_tests")
if(AKANTU_CONTACT)
  add_akantu_test(test_contact "test_contact")
endif(AKANTU_CONTACT)

#===============================================================================
add_mesh(test_solid_mechanics_model_square_mesh square.geo 2 1)
add_mesh(test_solid_mechanics_model_circle_mesh1 circle.geo 2 1 OUTPUT circle1.msh)
add_mesh(test_solid_mechanics_model_circle_mesh2 circle.geo 2 2 OUTPUT circle2.msh)

register_test(test_solid_mechanics_model_square test_solid_mechanics_model_square.cc)
add_dependencies(test_solid_mechanics_model_square
  test_solid_mechanics_model_square_mesh
  )

register_test(test_solid_mechanics_model_circle_2 test_solid_mechanics_model_circle_2.cc)
add_dependencies(test_solid_mechanics_model_circle_2
  test_solid_mechanics_model_circle_mesh2)

#===============================================================================
add_mesh(test_bar_traction_2d_mesh1 bar.geo 2 1 OUTPUT bar1.msh)
add_mesh(test_bar_traction_2d_mesh2 bar.geo 2 2 OUTPUT bar2.msh)
add_mesh(test_bar_traction_2d_mesh_structured1 bar_structured.geo 2 1 OUTPUT bar_structured1.msh)

register_test(test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d.cc)
add_dependencies(test_solid_mechanics_model_bar_traction2d test_bar_traction_2d_mesh1 test_bar_traction_2d_mesh2)

register_test(test_solid_mechanics_model_bar_traction2d_structured
  test_solid_mechanics_model_bar_traction2d_structured.cc)
add_dependencies(test_solid_mechanics_model_bar_traction2d_structured test_bar_traction_2d_mesh_structured1)

if(AKANTU_SCOTCH AND AKANTU_MPI)
  register_test(test_solid_mechanics_model_bar_traction2d_parallel
    test_solid_mechanics_model_bar_traction2d_parallel.cc)
  add_dependencies(test_solid_mechanics_model_bar_traction2d_parallel
    test_bar_traction_2d_mesh2)

  add_mesh(test_solid_mechanics_model_segment_mesh segment.geo 1 2)

  register_test(test_solid_mechanics_model_segment_parallel
    test_solid_mechanics_model_segment_parallel.cc)
  add_dependencies(test_solid_mechanics_model_segment_parallel
    test_solid_mechanics_model_segment_mesh)
endif()

if(AKANTU_MUMPS)
  register_test(test_solid_mechanics_model_bar_traction2d_mass_not_lumped
    test_solid_mechanics_model_bar_traction2d_mass_not_lumped.cc)
  add_dependencies(test_solid_mechanics_model_bar_traction2d_mass_not_lumped test_bar_traction_2d_mesh1 test_bar_traction_2d_mesh2)

  add_mesh(test_implicit_mesh1 square_implicit.geo 2 1 OUTPUT square_implicit1.msh)
  add_mesh(test_implicit_mesh2 square_implicit.geo 2 2 OUTPUT square_implicit2.msh)
  register_test(test_solid_mechanics_model_implicit_2d test_solid_mechanics_model_implicit_2d.cc)
  add_dependencies(test_solid_mechanics_model_implicit_2d test_implicit_mesh1 test_implicit_mesh2)

  add_mesh(test_solid_mechanics_model_segment_mesh1 segment.geo 1 1 OUTPUT segment1.msh)
  register_test(test_solid_mechanics_model_implicit_1d test_solid_mechanics_model_implicit_1d.cc)
  add_dependencies(test_solid_mechanics_model_implicit_1d test_solid_mechanics_model_segment_mesh1)

  add_mesh(test_implicit_beam_2d_1 beam_2d.geo 2 1 OUTPUT beam_2d_lin.msh)
  add_mesh(test_implicit_beam_2d_2 beam_2d.geo 2 2 OUTPUT beam_2d_quad.msh)
  add_mesh(test_implicit_beam_3d_2 beam_3d.geo 3 2 OUTPUT beam_3d_quad.msh)
  add_mesh(test_implicit_beam_3d_1 beam_3d.geo 3 1 OUTPUT beam_3d_lin.msh)
  register_test(test_solid_mechanics_model_implicit_dynamic_2d
    test_solid_mechanics_model_implicit_dynamic_2d.cc)
  add_dependencies(test_solid_mechanics_model_implicit_dynamic_2d test_implicit_beam_2d_1 test_implicit_beam_2d_2 test_implicit_beam_3d_2 test_implicit_beam_3d_1)
  copy_files(test_solid_mechanics_model_implicit_dynamic_2d material_implicit_dynamic.dat)
endif()

#===============================================================================
add_mesh(test_cube3d_mesh1 cube.geo 3 1 OUTPUT cube1.msh)
add_mesh(test_cube3d_mesh2 cube.geo 3 2 OUTPUT cube2.msh)

register_test(test_solid_mechanics_model_cube3d
  test_solid_mechanics_model_cube3d.cc)
add_dependencies(test_solid_mechanics_model_cube3d test_cube3d_mesh1)

register_test(test_solid_mechanics_model_cube3d_tetra10
  test_solid_mechanics_model_cube3d_tetra10.cc)
add_dependencies(test_solid_mechanics_model_cube3d_tetra10 test_cube3d_mesh2)

#===============================================================================
file(COPY material.dat DESTINATION .)
file(COPY test_cst_energy.pl DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
