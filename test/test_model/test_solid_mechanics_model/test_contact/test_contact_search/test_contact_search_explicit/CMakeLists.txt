#===============================================================================
# @file   CMakeLists.txt
#
# @author David Simon Kammer <david.kammer@epfl.ch>
#
# @date   Mon Jan 17 11:53:13 2011
#
# @brief  configuration for explicit ContactSearch tests
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

set(LIST_TYPES
  triangle_3
  quadrangle_4
  tetrahedron_4
  #hexahedron_8
  )

#===============================================================================
foreach(_type ${LIST_TYPES})
  if(AKANTU_IOHELPER)
    register_test(test_contact_search_explicit_${_type} test_contact_search_explicit.cc io_helper_tools.cc)
  else()
    register_test(test_contact_search_explicit_${_type} test_contact_search_explicit.cc)
  endif()
  set_target_properties(test_contact_search_explicit_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ./data/${_type}.msh DESTINATION .)
  file(COPY ./ref_files/test_contact_search_explicit_${_type}.ref DESTINATION .)
endforeach()

#===============================================================================
file(COPY material.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
