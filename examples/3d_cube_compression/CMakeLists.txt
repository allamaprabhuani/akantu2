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
add_mesh(3d_cube_compression_mesh cube.geo 3 1)

register_example(3d_cube_compression
  3d_cube_compression.cc)
add_dependencies(3d_cube_compression
  3d_cube_compression_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)