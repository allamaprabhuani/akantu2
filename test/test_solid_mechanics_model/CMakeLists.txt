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

REGISTER_TEST(test_solid_mechanics_model_bar_traction2d test_solid_mechanics_model_bar_traction2d.cc)
add_mesh(bar bar.geo 2 1)
add_dependencies(test_solid_mechanics_model_bar_traction2d bar)

REGISTER_TEST(test_solid_mechanics_model test_solid_mechanics_model.cc)
add_mesh(triangle triangle.geo 2 1)
add_dependencies(test_solid_mechanics_model triangle)

REGISTER_TEST(test_solid_mechanics_model_cube3d test_solid_mechanics_model_cube3d.cc)
add_mesh(lcube cube.geo 3 1)
add_dependencies(test_solid_mechanics_model_cube3d lcube)

