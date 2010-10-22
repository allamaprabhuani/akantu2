#===============================================================================
# @file   CMakeLists.txt
# @author Leo 
# @date   Fri Sep 29 16:46:30 2010 ...più precisi di così
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(2d_compression_mesh square_2d.geo 2 1)

register_test(2d_compression
  2d_compression.cc)
add_dependencies(2d_compression 
  2d_compression_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )