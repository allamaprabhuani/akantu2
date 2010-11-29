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
add_mesh(test_facet_square_mesh square.geo 2 1)

register_test(test_facet_extraction_triangle_3 test_facet_extraction_triangle_3.cc)
add_dependencies(test_facet_extraction_triangle_3 test_facet_square_mesh)

#===============================================================================
add_mesh(test_facet_cube_mesh cube.geo 3 1)

register_test(test_facet_extraction_tetrahedron_4 test_facet_extraction_tetrahedron_4.cc)
add_dependencies(test_facet_extraction_tetrahedron_4 test_facet_cube_mesh)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)