#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Aug  8 15:39:24 2012
#
# @brief  
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
#===============================================================================

add_mesh(test_material_standard_linear_solid_deviatoric_relaxation_mesh test_material_standard_linear_solid_deviatoric_relaxation.geo 2 1)
register_test(test_material_standard_linear_solid_deviatoric_relaxation test_material_standard_linear_solid_deviatoric_relaxation.cc)
add_dependencies(test_material_standard_linear_solid_deviatoric_relaxation test_material_standard_linear_solid_deviatoric_relaxation_mesh)

register_test(test_material_standard_linear_solid_deviatoric_relaxation_tension test_material_standard_linear_solid_deviatoric_relaxation_tension.cc)
add_dependencies(test_material_standard_linear_solid_deviatoric_relaxation_tension test_material_standard_linear_solid_deviatoric_relaxation_mesh)
file(COPY material_standard_linear_solid_deviatoric_relaxation.dat DESTINATION .)

#===============================================================================
#foreach(_type ${LIST_TYPES})
#  if(AKANTU_IOHELPER)
#    register_test(test_material_stiffness_proportional_${_type} test_material_stiffness_proportional.cc io_helper_tools.cc)
#  else()
#    register_test(test_material_stiffness_proportional_${_type} test_material_stiffness_proportional.cc)
#  endif()
#  set_target_properties(test_material_stiffness_proportional_${_type}
#    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
#  file(COPY ./data/test_material_stiffness_proportional_${_type}.msh DESTINATION .)
#  file(COPY ./ref_files/test_material_stiffness_proportional_${_type}.ref DESTINATION .)
#endforeach()
#
##===============================================================================
#foreach(_type ${LIST_TYPES})
#  if(AKANTU_IOHELPER)
#    register_test(test_material_stiffness_proportional_damping_${_type} test_material_stiffness_proportional_damping.cc io_helper_tools.cc)
#  else()
#    register_test(test_material_stiffness_proportional_damping_${_type} test_material_stiffness_proportional_damping.cc)
#  endif()
#  set_target_properties(test_material_stiffness_proportional_damping_${_type}
#    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
#  file(COPY ./data/single_${_type}.msh DESTINATION .)
#endforeach()
#
##===============================================================================
#
#file(COPY material_stiffness_proportional.dat DESTINATION .)
#file(COPY material_stiffness_proportional_damping.dat DESTINATION .)
