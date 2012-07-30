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
macro(register_patch_test name type plane_strain)
  register_test(patch_test_explicit_${name} patch_test_explicit.cc)
  file(COPY ../data/${type}.msh DESTINATION .)
  set_target_properties(patch_test_explicit_${name}
    PROPERTIES COMPILE_DEFINITIONS "TYPE=_${type};PLANE_STRAIN=${plane_strain}")
endmacro()

foreach(_type ${LIST_TYPES})
  register_patch_test(${_type} ${_type} true)
endforeach()

foreach(_type ${LIST_TYPES_2D})
  register_patch_test(plane_strain_${_type} ${_type} true)
  register_patch_test(plane_stress_${_type} ${_type} false)
endforeach()

#===============================================================================
file(COPY ../data/material_check_stress_plane_strain.dat DESTINATION .)
file(COPY ../data/material_check_stress_plane_stress.dat DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)