#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date   Fri Sep 03 15:56:15 2010
#
# @brief  configuration for FEM tests
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
macro(register_fem_test operation type)
  set(_target test_${operation}_${_type})
  register_test(${_target}
    SOURCES test_${operation}.cc
    FILES_TO_COPY _${type}.msh
    COMPILE_OPTIONS TYPE=_${_type}
    PACKAGE core
    )
endmacro()
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

foreach(_type ${LIST_TYPES})
  register_fem_test(interpolate ${_type})
  register_fem_test(gradient    ${_type})
  register_fem_test(integrate   ${_type})
  register_fem_test(inverse_map ${_type})
endforeach()

#register_test(test_interpolate_bernoulli_beam_2 test_interpolate_bernoulli_beam_2.cc)
#add_mesh(test_fem_circle_1_mesh circle.geo 2 1 OUTPUT circle1.msh)
#add_mesh(test_fem_circle_2_mesh circle.geo 2 2 OUTPUT circle2.msh)
