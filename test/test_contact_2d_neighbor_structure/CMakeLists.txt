#===============================================================================
# @file   CMakeLists.txt
# @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
# @date   Thu Dec  9 10:19:20 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_contact_2d_neighbor_structure_mesh squares.geo 2 1)

register_test(test_contact_2d_neighbor_structure test_contact_2d_neighbor_structure.cc)
add_dependencies(test_contact_2d_neighbor_structure test_contact_2d_neighbor_structure_mesh)
