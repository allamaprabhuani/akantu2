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

add_subdirectory("test_materials")
#===============================================================================
add_mesh(test_solid_mechanics_model_square_mesh triangle.geo 2 1)
add_mesh(test_solid_mechanics_model_circle_mesh circle.geo 2 2)

register_test(test_solid_mechanics_model_square test_solid_mechanics_model_square.cc)
add_dependencies(test_solid_mechanics_model_square
  test_solid_mechanics_model_square_mesh
  test_solid_mechanics_model_circle_mesh)

register_test(test_solid_mechanics_model_circle_2 test_solid_mechanics_model_circle_2.cc)

#===============================================================================
add_mesh(test_bar_traction_2d_mesh1 bar.geo 2 1 OUTPUT bar1.msh)
add_mesh(test_bar_traction_2d_mesh2 bar.geo 2 2 OUTPUT bar2.msh)

register_test(test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d.cc)
add_dependencies(test_solid_mechanics_model_bar_traction2d test_bar_traction_2d_mesh1 test_bar_traction_2d_mesh2)

if(AKANTU_SCOTCH_ON AND AKANTU_MPI_ON)
  register_test(test_solid_mechanics_model_bar_traction2d_parallel
    test_solid_mechanics_model_bar_traction2d_parallel.cc)
  add_dependencies(test_solid_mechanics_model_bar_traction2d_parallel
    test_bar_traction_2d_mesh2)


  add_mesh(test_solid_mechanics_model_segment_mesh segment.geo 1 2)

  register_test(test_solid_mechanics_model_segment_parallel
    test_solid_mechanics_model_segment_parallel.cc)
  add_dependencies(test_solid_mechanics_model_segment_parallel
    test_solid_mechanics_model_segment_mesh)
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

