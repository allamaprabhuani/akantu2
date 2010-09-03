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

REGISTER_TEST(test_interpolate_line_1 test_interpolate_line_1.cc)
add_mesh(line line.geo 1 1)
add_dependencies(test_interpolate_line_1 line)

REGISTER_TEST(test_interpolate_triangle_1 test_interpolate_triangle_1.cc)
add_mesh(ltriangle triangle.geo 2 1)
add_dependencies(test_interpolate_triangle_1 ltriangle)

REGISTER_TEST(test_interpolate_tetrahedra_1 test_interpolate_tetrahedra_1.cc)
add_mesh(llcube cube.geo 3 1)
add_dependencies(test_interpolate_tetrahedra_1 llcube)

