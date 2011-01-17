#===============================================================================
# @file   CMakeLists.txt
# @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
# @date   Fri Dec  3 15:12:22 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_contact_2d_expli_mesh squares.geo 2 1)

register_test(test_contact_2d_expli test_contact_2d_expli.cc)
add_dependencies(test_contact_2d_expli test_contact_2d_expli_mesh)
