#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Tue Oct 26 16:40:24 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_regular_grid_2d_mesh squares.geo 2 1)

register_test(test_contact_regular_grid_triangle_3 test_contact_regular_grid_triangle_3.cc)
add_dependencies(test_contact_regular_grid_triangle_3 test_regular_grid_2d_mesh)

#===============================================================================
add_mesh(test_regular_grid_3d_mesh cubes.geo 3 1)

register_test(test_contact_regular_grid test_contact_regular_grid.cc)
add_dependencies(test_contact_regular_grid test_regular_grid_3d_mesh)

#===============================================================================
#add_mesh(test_regular_grid_3d_mesh cubes.geo 3 1)

register_test(test_contact_regular_grid_nodes test_contact_regular_grid_nodes.cc)
add_dependencies(test_contact_regular_grid_nodes test_regular_grid_3d_mesh)