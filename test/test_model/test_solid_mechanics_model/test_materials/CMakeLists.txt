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

add_mesh(test_local_material_square_mesh triangle.geo 2 2)
add_mesh(test_local_material_circle_mesh circle.geo 2 2)
add_mesh(test_local_material_bar_mesh bar.geo 2 2)
add_mesh(test_local_material_barre_trou_mesh barre_trou.geo 2 2)

register_test(test_local_material material_damage.cc test_local_material.cc)
add_dependencies(test_local_material
  test_local_material_square_mesh
  test_local_material_circle_mesh
  test_local_material_bar_mesh
  test_local_material_barre_trou_mesh
  )

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)