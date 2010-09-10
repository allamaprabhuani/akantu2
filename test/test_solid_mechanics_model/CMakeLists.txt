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

register_test(test_solid_mechanics_model_bar_traction2d test_solid_mechanics_model_bar_traction2d.cc)
add_mesh(bar bar.geo 2 1)
add_dependencies(test_solid_mechanics_model_bar_traction2d bar)

register_test(test_solid_mechanics_model test_solid_mechanics_model.cc)
add_mesh(triangle triangle.geo 2 1)
add_dependencies(test_solid_mechanics_model triangle)

register_test(test_solid_mechanics_model_cube3d test_solid_mechanics_model_cube3d.cc)
add_mesh(cube_mesh cube.geo 3 1)
add_dependencies(test_solid_mechanics_model_cube3d cube_mesh)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )