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

REGISTER_TEST(test_mesh_io_msh test_mesh_io_msh.cc)
add_mesh(test_mesh_io_mesh cube.geo 3 1)
add_dependencies(test_mesh_io_msh test_mesh_io_mesh)