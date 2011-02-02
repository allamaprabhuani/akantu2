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

add_mesh(test_solver_mesh triangle.geo 2 1)

register_test(test_sparse_matrix_profile test_sparse_matrix_profile.cc)
add_dependencies(test_sparse_matrix_profile test_solver_mesh)

register_test(test_sparse_matrix_assemble test_sparse_matrix_assemble.cc)
add_dependencies(test_sparse_matrix_assemble test_solver_mesh)

register_test(test_solver_mumps test_solver_mumps.cc)