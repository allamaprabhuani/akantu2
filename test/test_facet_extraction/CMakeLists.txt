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

REGISTER_TEST(test_facet_extraction_tetra1 test_facet_extraction_tetra1.cc)
add_mesh(cube cube.geo 3 1)
add_dependencies(test_facet_extraction_tetra1 cube)

REGISTER_TEST(test_facet_extraction_triangle1 test_facet_extraction_triangle1.cc)
add_mesh(square square.geo 2 1)
add_dependencies(test_facet_extraction_triangle1 square)

