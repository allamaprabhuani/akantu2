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
add_mesh(test_fem_segment_2_mesh segment.geo 1 1 OUTPUT segment_2.msh)
add_mesh(test_fem_segment_3_mesh segment.geo 1 2 OUTPUT segment_3.msh)

#===============================================================================
# 1st order
register_test(test_integrate_segment_2 test_integrate_segment_2.cc)
add_dependencies(test_integrate_segment_2 test_fem_segment_2_mesh)
#===============================================================================
register_test(test_interpolate_segment_2 test_interpolate_segment_2.cc)
add_dependencies(test_interpolate_segment_2 test_fem_segment_2_mesh)
#===============================================================================
register_test(test_gradient_segment_2 test_gradient_segment_2.cc)
add_dependencies(test_gradient_segment_2 test_fem_segment_2_mesh)

#===============================================================================
# 2nd order
register_test(test_integrate_segment_3 test_integrate_segment_3.cc)
add_dependencies(test_integrate_segment_3 test_fem_segment_3_mesh)
#===============================================================================
register_test(test_interpolate_segment_3 test_interpolate_segment_3.cc)
add_dependencies(test_interpolate_segment_3 test_fem_segment_3_mesh)
#===============================================================================
register_test(test_gradient_segment_3 test_gradient_segment_3.cc)
add_dependencies(test_gradient_segment_3 test_fem_segment_3_mesh)

#===============================================================================
# 2D
#===============================================================================
add_mesh(test_fem_triangle_3_mesh triangle.geo 2 1 OUTPUT triangle_3.msh)
add_mesh(test_fem_triangle_6_mesh triangle.geo 2 2 OUTPUT triangle_6.msh)

#===============================================================================
# 1st order
register_test(test_integrate_triangle_3 test_integrate_triangle_3.cc)
add_dependencies(test_integrate_triangle_3 test_fem_triangle_3_mesh)
#===============================================================================
register_test(test_interpolate_triangle_3 test_interpolate_triangle_3.cc)
add_dependencies(test_interpolate_triangle_3 test_fem_triangle_3_mesh)
#===============================================================================
register_test(test_gradient_triangle_3 test_gradient_triangle_3.cc)
add_dependencies(test_gradient_triangle_3 test_fem_triangle_3_mesh)

#===============================================================================
# 2nd order
register_test(test_integrate_triangle_6 test_integrate_triangle_6.cc)
add_dependencies(test_integrate_triangle_6 test_fem_triangle_6_mesh)
#===============================================================================
register_test(test_interpolate_triangle_6 test_interpolate_triangle_6.cc)
add_dependencies(test_interpolate_triangle_6 test_fem_triangle_6_mesh)
#===============================================================================
register_test(test_gradient_triangle_6 test_gradient_triangle_6.cc)
add_dependencies(test_gradient_triangle_6 test_fem_triangle_6_mesh)

#===============================================================================
# 3D
#===============================================================================
add_mesh(test_fem_tetrahedron_6_mesh cube.geo 3 1)

register_test(test_integrate_tetrahedron_6 test_integrate_tetrahedron_6.cc)
add_dependencies(test_integrate_tetrahedron_6 test_fem_tetrahedron_6_mesh)
#===============================================================================
register_test(test_interpolate_tetrahedron_6 test_interpolate_tetrahedron_6.cc)
add_dependencies(test_interpolate_tetrahedron_6 test_fem_tetrahedron_6_mesh)
#===============================================================================
register_test(test_gradient_tetrahedron_6 test_gradient_tetrahedron_6.cc)
add_dependencies(test_gradient_tetrahedron_6 test_fem_tetrahedron_6_mesh)


