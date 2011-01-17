#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Fri Dec 03 12:02:24 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
#add_mesh(test_regular_grid_2d_mesh squares.geo 2 1)

#register_test(test_regular_grid_triangle_3 test_regular_grid_triangle_3.cc)
#add_dependencies(test_regular_grid_triangle_3 test_regular_grid_2d_mesh)

#===============================================================================
add_mesh(test_3d_explicit_mesh cubes.geo 3 1)

register_test(test_3d_explicit_tetrahedron_4 test_3d_explicit_tetrahedron_4.cc)
add_dependencies(test_3d_explicit_tetrahedron_4 test_3d_explicit_mesh)

#===============================================================================
#add_mesh(test_regular_grid_3d_mesh cubes.geo 3 1)

#register_test(test_regular_grid_tetrahedron_4_nodes test_regular_grid_tetrahedron_4_nodes.cc)
#add_dependencies(test_regular_grid_tetrahedron_4_nodes test_regular_grid_3d_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
