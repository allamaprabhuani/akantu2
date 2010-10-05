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

add_mesh(test_mesh_partitionate_mesh
  triangle.geo 2 2)

register_test(test_mesh_partitionate_scotch
  test_mesh_partitionate_scotch.cc)
add_dependencies(test_mesh_partitionate_scotch
  test_mesh_partitionate_mesh)