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
# 1D
#===============================================================================
add_mesh(test_fem_line_1_mesh line.geo 1 1 OUTPUT line1.msh)
add_mesh(test_fem_line_2_mesh line.geo 1 2 OUTPUT line2.msh)

#===============================================================================
# 1st order
register_test(test_integrate_segment_2 test_integrate_segment_2.cc)
add_dependencies(test_integrate_segment_2 test_fem_line_1_mesh)
#===============================================================================
register_test(test_interpolate_segment_2 test_interpolate_segment_2.cc)
add_dependencies(test_interpolate_segment_2 test_fem_line_1_mesh)
#===============================================================================
register_test(test_gradient_segment_2 test_gradient_segment_2.cc)
add_dependencies(test_gradient_segment_2 test_fem_line_1_mesh)

#===============================================================================
# 2nd order
register_test(test_integrate_segment_3 test_integrate_segment_3.cc)
add_dependencies(test_integrate_segment_3 test_fem_line_2_mesh)
#===============================================================================
register_test(test_interpolate_segment_3 test_interpolate_segment_3.cc)
add_dependencies(test_interpolate_segment_3 test_fem_line_2_mesh)
#===============================================================================
register_test(test_gradient_segment_3 test_gradient_segment_3.cc)
add_dependencies(test_gradient_segment_3 test_fem_line_2_mesh)

#===============================================================================
# 2D
#===============================================================================
add_mesh(test_fem_square_1_mesh square.geo 2 1 OUTPUT square1.msh)
add_mesh(test_fem_square_2_mesh square.geo 2 2 OUTPUT square2.msh)
add_mesh(test_fem_circle_1_mesh circle.geo 2 1 OUTPUT circle1.msh)
add_mesh(test_fem_circle_2_mesh circle.geo 2 2 OUTPUT circle2.msh)
add_mesh(test_fem_square_struct_1_mesh square_structured.geo 2 1 OUTPUT square_structured1.msh)

#===============================================================================
# 1st order unstructured
register_test(test_integrate_triangle_3 test_integrate_triangle_3.cc)
add_dependencies(test_integrate_triangle_3
  test_fem_square_1_mesh
  test_fem_circle_1_mesh)
#===============================================================================
register_test(test_interpolate_triangle_3 test_interpolate_triangle_3.cc)
add_dependencies(test_interpolate_triangle_3 test_fem_square_1_mesh)
#===============================================================================
register_test(test_gradient_triangle_3 test_gradient_triangle_3.cc)
add_dependencies(test_gradient_triangle_3 test_fem_square_1_mesh)
#===============================================================================
# 1st order structured
register_test(test_integrate_quadrangle_4 test_integrate_quadrangle_4.cc)
add_dependencies(test_integrate_quadrangle_4
  test_fem_square_1_mesh
  test_fem_circle_1_mesh)
#===============================================================================
register_test(test_interpolate_quadrangle_4 test_interpolate_quadrangle_4.cc)
add_dependencies(test_interpolate_quadrangle_4 test_fem_square_struct_1_mesh)
#===============================================================================
register_test(test_gradient_quadrangle_4 test_gradient_quadrangle_4.cc)
add_dependencies(test_gradient_quadrangle_4 test_fem_square_struct_1_mesh)


#===============================================================================
# 2nd order
register_test(test_integrate_triangle_6 test_integrate_triangle_6.cc)
add_dependencies(test_integrate_triangle_6
  test_fem_square_2_mesh
  test_fem_circle_2_mesh)
#===============================================================================
register_test(test_interpolate_triangle_6 test_interpolate_triangle_6.cc)
add_dependencies(test_interpolate_triangle_6 test_fem_square_2_mesh)
#===============================================================================
register_test(test_gradient_triangle_6 test_gradient_triangle_6.cc)
add_dependencies(test_gradient_triangle_6 test_fem_square_2_mesh)

#===============================================================================
# 3D
#===============================================================================
add_mesh(test_fem_cube_1_mesh cube.geo 3 1 OUTPUT cube1.msh)
add_mesh(test_fem_cube_2_mesh cube.geo 3 2 OUTPUT cube2.msh)

#===============================================================================
# 1st order
register_test(test_integrate_tetrahedron_4 test_integrate_tetrahedron_4.cc)
add_dependencies(test_integrate_tetrahedron_4 test_fem_cube_1_mesh)
#===============================================================================
register_test(test_interpolate_tetrahedron_4 test_interpolate_tetrahedron_4.cc)
add_dependencies(test_interpolate_tetrahedron_4 test_fem_cube_1_mesh)
#===============================================================================
register_test(test_gradient_tetrahedron_4 test_gradient_tetrahedron_4.cc)
add_dependencies(test_gradient_tetrahedron_4 test_fem_cube_1_mesh)

#===============================================================================
# 2nd order
register_test(test_integrate_tetrahedron_10 test_integrate_tetrahedron_10.cc)
add_dependencies(test_integrate_tetrahedron_10 test_fem_cube_2_mesh)
#===============================================================================
register_test(test_interpolate_tetrahedron_10 test_interpolate_tetrahedron_10.cc)
add_dependencies(test_interpolate_tetrahedron_10 test_fem_cube_2_mesh)
#===============================================================================
register_test(test_gradient_tetrahedron_10 test_gradient_tetrahedron_10.cc)
add_dependencies(test_gradient_tetrahedron_10 test_fem_cube_2_mesh)
