#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(test_synchronizer_communication_mesh
  triangle.geo 2 1)

register_test(test_synchronizer_communication
  test_synchronizer_communication.cc)
add_dependencies(test_synchronizer_communication
  test_synchronizer_communication_mesh)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/test_synchronizer_communication.sh
  ${CMAKE_CURRENT_BINARY_DIR}/test_synchronizer_communication.sh
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)