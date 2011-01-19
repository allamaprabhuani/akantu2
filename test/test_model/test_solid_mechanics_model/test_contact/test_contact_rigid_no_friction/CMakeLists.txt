#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Wed Jan 19 12:37:24 2011
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
add_mesh(test_squares_mesh squares.geo 2 1)

register_test(test_contact_rigid_no_friction_triangle_3 test_contact_rigid_no_friction_triangle_3.cc)
add_dependencies(test_contact_rigid_no_friction_triangle_3 test_squares_mesh)

#===============================================================================
add_mesh(test_cubes_mesh cubes.geo 3 1)

register_test(test_contact_rigid_no_friction_tetrahedron_4 test_contact_rigid_no_friction_tetrahedron_4.cc)
add_dependencies(test_contact_rigid_no_friction_tetrahedron_4 test_cubes_mesh)

#===============================================================================
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/material.dat
  ${CMAKE_CURRENT_BINARY_DIR}/material.dat
  COPYONLY
  )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
