#===============================================================================
# @file   rve_tools.cmake
#
# @author Emil Gallyamov <emil.gallyamov@epfl.ch>
#
# @date creation: Mon Feb 22 2022
#
# @brief  package description for rve tools
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

package_declare(rve_tools DEFAULT OFF
  DESCRIPTION "Multi-scale model package"
  DEPENDS extra_materials cohesive_element)

package_declare_sources(rve_tools
  model/rve_tools/rve_tools.cc
  model/rve_tools/rve_tools.hh
  model/rve_tools/materials/material_FE2/material_FE2.hh
  model/rve_tools/materials/material_FE2/material_FE2.cc
  model/rve_tools/solid_mechanics_model_RVE.hh
  model/rve_tools/solid_mechanics_model_RVE.cc
  model/rve_tools/mesh_utils/nodes_flag_updater.hh
  model/rve_tools/mesh_utils/nodes_flag_updater.cc
  model/rve_tools/mesh_utils/nodes_flag_updater_inline_impl.hh
  model/rve_tools/mesh_utils/nodes_eff_stress_updater.hh
  model/rve_tools/mesh_utils/nodes_eff_stress_updater.cc
  model/rve_tools/mesh_utils/nodes_eff_stress_updater_inline_impl.hh
  model/rve_tools/mesh_utils/crack_numbers_updater.hh
  model/rve_tools/mesh_utils/crack_numbers_updater.cc
  model/rve_tools/mesh_utils/crack_numbers_updater_inline_impl.hh
  model/rve_tools/materials/material_cohesive_linear_sequential.hh
  model/rve_tools/materials/material_cohesive_linear_sequential.cc
  model/rve_tools/materials/material_cohesive_linear_sequential_inline_impl.hh
)
