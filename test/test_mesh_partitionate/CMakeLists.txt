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

register_test(test_mesh_partitionate_scotch
  test_mesh_partitionate.cc)
add_mesh(test_mesh_partitionate_mesh
  triangle.geo 2 1)
add_dependencies(test_mesh_partitionate_scotch
  test_mesh_partitionate_mesh)