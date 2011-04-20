#===============================================================================
# @file   CMakeLists.txt
# @author David Kammer <david.kammer@epfl.ch>
# @date   Wed Feb 16 09:46:59 2010
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
  segment_2
  segment_3
  triangle_3
  triangle_6
  quadrangle_4
  tetrahedron_4
  tetrahedron_10
  hexahedron_8
  )

foreach(_type ${LIST_TYPES})
  register_test(test_check_stress_${_type} test_check_stress.cc)
  set_target_properties(test_check_stress_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ${_type}.msh DESTINATION .)
endforeach()

#===============================================================================
file(COPY material_check_stress.dat DESTINATION .)


