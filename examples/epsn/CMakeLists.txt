#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Dec 8 2010 16:46:30 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(epsn_bar_mesh bar.geo 2 1)

register_example(akantu_epsn_example
  akantu_epsn_example.cc)
add_dependencies(akantu_epsn_example
  epsn_bar_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/bar_traction.xml
  ${CMAKE_CURRENT_BINARY_DIR}/bar_traction.xml
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)