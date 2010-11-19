#===============================================================================
# @file   CMakeLists.txt
# @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
# @date   Mon Oct 25 09:40:24 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_surface_extraction_2d_mesh squares.geo 2 1)

register_test(test_surface_extraction_triangle_3 test_surface_extraction_triangle_3.cc)
add_dependencies(test_surface_extraction_triangle_3 test_surface_extraction_2d_mesh)

#===============================================================================
add_mesh(test_surface_extraction_3d_mesh cubes.geo 3 1)

register_test(test_surface_extraction_tetrahedron_6 test_surface_extraction_tetrahedron_6.cc)
add_dependencies(test_surface_extraction_tetrahedron_6 test_surface_extraction_3d_mesh)