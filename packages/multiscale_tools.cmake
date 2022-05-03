#===============================================================================
# @file   multiscale_tools.cmake
#
# @author Emil Gallyamov <emil.gallyamov@epfl.ch>
#
# @date creation: Mon Feb 22 2022
#
# @brief  package description for multiscale tools
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

package_declare(multiscale_tools
  DESCRIPTION "Multi-scale model package"
  DEPENDS extra_materials cohesive_element)

package_declare_sources(multiscale_tools
  model/multiscale_tools/rve_tools.cc
  model/multiscale_tools/rve_tools.hh
  model/multiscale_tools/expanding_material_selector.hh
  model/multiscale_tools/bc_functors.hh
  model/multiscale_tools/materials/material_FE2/material_FE2.hh
  model/multiscale_tools/materials/material_FE2/material_FE2.cc
  model/multiscale_tools/solid_mechanics_model_RVE.hh
  model/multiscale_tools/solid_mechanics_model_RVE.cc
  model/multiscale_tools/mesh_utils/nodes_flag_updater.hh
  model/multiscale_tools/mesh_utils/nodes_flag_updater.cc
  model/multiscale_tools/mesh_utils/nodes_flag_updater_inline_impl.hh
  model/multiscale_tools/mesh_utils/nodes_eff_stress_updater.hh
  model/multiscale_tools/mesh_utils/nodes_eff_stress_updater.cc
  model/multiscale_tools/mesh_utils/nodes_eff_stress_updater_inline_impl.hh
  model/multiscale_tools/mesh_utils/crack_numbers_updater.hh
  model/multiscale_tools/mesh_utils/crack_numbers_updater.cc
  model/multiscale_tools/mesh_utils/crack_numbers_updater_inline_impl.hh
  model/multiscale_tools/materials/material_cohesive_linear_sequential.hh
  model/multiscale_tools/materials/material_cohesive_linear_sequential.cc
  model/multiscale_tools/materials/material_cohesive_linear_sequential_inline_impl.hh
)
