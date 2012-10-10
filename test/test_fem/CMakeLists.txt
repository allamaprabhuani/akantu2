#===============================================================================
# @file   CMakeLists.txt
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
#register_test(test_interpolate_bernoulli_beam_2 test_interpolate_bernoulli_beam_2.cc)

#===============================================================================
set(LIST_TYPES
  segment_2
  segment_3
  triangle_3
  triangle_6
  quadrangle_4
  quadrangle_8
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  )

add_mesh(test_fem_segment_2_mesh      line.geo              1 1 OUTPUT _segment_2.msh)
add_mesh(test_fem_segment_3_mesh      line.geo              1 2 OUTPUT _segment_3.msh)
add_mesh(test_fem_triangle_3_mesh     square.geo            2 1 OUTPUT _triangle_3.msh)
add_mesh(test_fem_triangle_6_mesh     square.geo            2 2 OUTPUT _triangle_6.msh)
add_mesh(test_fem_quadrangle_4_mesh   square_structured.geo 2 1 OUTPUT _quadrangle_4.msh)
add_mesh(test_fem_quadrangle_8_mesh   square_structured.geo 2 2 OUTPUT _quadrangle_8.msh)
add_mesh(test_fem_tetrahedron_4_mesh  cube.geo              3 1 OUTPUT _tetrahedron_4.msh)
add_mesh(test_fem_tetrahedron_10_mesh cube.geo              3 2 OUTPUT _tetrahedron_10.msh)
add_mesh(test_fem_hexahedron_8_mesh   hexa_structured.geo   3 1 OUTPUT _hexahedron_8.msh)


macro(register_fem_test operation type)
  set(_target test_${operation}_${_type})
  register_test(${_target} test_${operation}.cc)
  set_target_properties(${_target}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  if(NOT EXISTS _${_type}.msh)
    add_dependencies(${_target} test_fem_${_type}_mesh)
  endif()
endmacro()


foreach(_type ${LIST_TYPES})
  register_fem_test(interpolate ${_type})
  register_fem_test(gradient    ${_type})
  register_fem_test(integrate   ${_type})
  register_fem_test(inverse_map ${_type})

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_type}.msh)
    file(COPY _${_type}.msh DESTINATION .)
  endif()
endforeach()

#add_mesh(test_fem_circle_1_mesh circle.geo 2 1 OUTPUT circle1.msh)
#add_mesh(test_fem_circle_2_mesh circle.geo 2 2 OUTPUT circle2.msh)
