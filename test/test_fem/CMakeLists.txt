#===============================================================================
# @file   CMakeLists.txt
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
# 1D
#===============================================================================
add_mesh(test_fem_line1_mesh line.geo 1 1)
register_test(test_integrate_line_1 test_integrate_line_1.cc)
add_dependencies(test_integrate_line_1 line)

#===============================================================================
register_test(test_interpolate_line_1 test_interpolate_line_1.cc)
add_dependencies(test_interpolate_line_1 test_fem_line1_mesh)
#===============================================================================
register_test(test_gradient_line_1 test_gradient_line_1.cc)
add_dependencies(test_gradient_line_1 test_fem_line1_mesh)
register_test(test_integrate_triangle_1 test_integrate_triangle_1.cc)
add_dependencies(test_integrate_triangle_1 ltriangle)


#===============================================================================
# 2D
#===============================================================================
add_mesh(test_fem_triangle1_mesh triangle.geo 2 1)
register_test(test_integrate_tetrahedra_1 test_integrate_tetrahedra_1.cc)
add_dependencies(test_integrate_tetrahedra_1 llcube)

#===============================================================================
register_test(test_interpolate_triangle_1 test_interpolate_triangle_1.cc)
add_dependencies(test_interpolate_triangle_1 test_fem_triangle1_mesh)
#===============================================================================
register_test(test_gradient_triangle_1 test_gradient_triangle_1.cc)
add_dependencies(test_gradient_triangle_1 test_fem_triangle1_mesh)

#===============================================================================
# 3D
#===============================================================================
add_mesh(test_fem_tetra1_mesh cube.geo 3 1)

#===============================================================================
register_test(test_interpolate_tetrahedra_1 test_interpolate_tetrahedra_1.cc)
add_dependencies(test_interpolate_tetrahedra_1 test_fem_tetra1_mesh)
#===============================================================================
register_test(test_gradient_tetrahedra_1 test_gradient_tetrahedra_1.cc)
add_dependencies(test_gradient_tetrahedra_1 test_fem_tetra1_mesh)


