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

add_mesh(square1 triangle.geo 2 1)
add_mesh(circle2 circle.geo 2 2)

register_test(test_local_material material_damage.cc test_local_material.cc) 
add_dependencies(test_local_material
  square1
  circle2)
#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

