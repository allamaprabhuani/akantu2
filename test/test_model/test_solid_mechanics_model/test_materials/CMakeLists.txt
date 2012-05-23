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

add_mesh(test_local_material_square_mesh triangle.geo 2 2)
add_mesh(test_local_material_circle_mesh circle.geo 2 2)
add_mesh(test_local_material_bar_mesh bar.geo 2 1)
add_mesh(test_local_material_barre_trou_mesh barre_trou.geo 2 2)
#add_mesh(test_square_structured_mesh square_structured.geo 2 1)
#add_mesh(test_single_triangle_mesh single_triangle.geo 2 1)

register_test(test_local_material local_material_damage.cc test_local_material.cc)
add_dependencies(test_local_material
  test_local_material_square_mesh
  test_local_material_circle_mesh

  test_local_material_barre_trou_mesh
  )

add_mesh(test_material_non_local_mesh mesh.geo 2 2 OUTPUT mesh.msh)
register_test(test_material_non_local test_material_non_local.cc)
add_dependencies(test_material_non_local
  test_material_non_local_mesh
  )

register_test(test_material_damage_non_local test_material_damage_non_local.cc)
add_dependencies(test_material_damage_non_local
  test_local_material_barre_trou_mesh
  )

#===============================================================================
file(COPY material.dat DESTINATION .)
file(COPY material_non_local.dat DESTINATION .)
file(COPY material_damage_non_local.dat DESTINATION .)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)

#===============================================================================
set(LIST_TYPES
  triangle_3
  quadrangle_4
  )

#===============================================================================
foreach(_type ${LIST_TYPES})
  if(AKANTU_IOHELPER)
    register_test(test_material_viscoelastic_${_type} test_material_viscoelastic.cc io_helper_tools.cc)
  else()
    register_test(test_material_viscoelastic_${_type} test_material_viscoelastic.cc)
  endif()
  set_target_properties(test_material_viscoelastic_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ./data/test_material_viscoelastic_${_type}.msh DESTINATION .)
  file(COPY ./ref_files/test_material_viscoelastic_${_type}.ref DESTINATION .)
endforeach()

add_mesh(test_material_viscoelastic_relaxation_mesh test_material_viscoelastic_relaxation.geo 2 1)
register_test(test_material_viscoelastic_relaxation test_material_viscoelastic_relaxation.cc)
add_dependencies(test_material_viscoelastic_relaxation test_material_viscoelastic_relaxation_mesh)
file(COPY material_viscoelastic_relaxation.dat DESTINATION .)

#===============================================================================
foreach(_type ${LIST_TYPES})
  if(AKANTU_IOHELPER)
    register_test(test_material_elastic_caughey_${_type} test_material_elastic_caughey.cc io_helper_tools.cc)
  else()
    register_test(test_material_elastic_caughey_${_type} test_material_elastic_caughey.cc)
  endif()
  set_target_properties(test_material_elastic_caughey_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ./data/test_material_elastic_caughey_${_type}.msh DESTINATION .)
  file(COPY ./ref_files/test_material_elastic_caughey_${_type}.ref DESTINATION .)
endforeach()

#===============================================================================
foreach(_type ${LIST_TYPES})
  if(AKANTU_IOHELPER)
    register_test(test_material_elastic_caughey_damping_${_type} test_material_elastic_caughey_damping.cc io_helper_tools.cc)
  else()
    register_test(test_material_elastic_caughey_damping_${_type} test_material_elastic_caughey_damping.cc)
  endif()
  set_target_properties(test_material_elastic_caughey_damping_${_type}
    PROPERTIES COMPILE_DEFINITIONS TYPE=_${_type})
  file(COPY ./data/single_${_type}.msh DESTINATION .)
endforeach()

#===============================================================================
file(COPY material_viscoelastic.dat DESTINATION .)
file(COPY material_elastic_caughey.dat DESTINATION .)
file(COPY material_elastic_caughey_damping.dat DESTINATION .)
