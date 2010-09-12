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
add_mesh(test_solid_mechanics_model_mesh triangle.geo 2 1)

register_test(test_solid_mechanics_model test_solid_mechanics_model.cc)
add_dependencies(test_solid_mechanics_model test_solid_mechanics_model_mesh)

#===============================================================================
add_mesh(test_bar_traction_2d_mesh bar.geo 2 1)

register_test(test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d.cc)
add_dependencies(test_solid_mechanics_model_bar_traction2d test_bar_traction_2d_mesh)

if(AKANTU_SCOTCH_ON AND AKANTU_MPI_ON)
  register_test(test_solid_mechanics_model_bar_traction2d_parallel
    test_solid_mechanics_model_bar_traction2d_parallel.cc)
  add_dependencies(test_solid_mechanics_model_bar_traction2d_parallel
    test_bar_traction_2d_mesh)
endif()

#===============================================================================
add_mesh(test_cube3d_mesh cube.geo 3 1)

register_test(test_solid_mechanics_model_cube3d
  test_solid_mechanics_model_cube3d.cc)
add_dependencies(test_solid_mechanics_model_cube3d test_cube3d_mesh)


#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )