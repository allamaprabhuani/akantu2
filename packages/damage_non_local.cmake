#===============================================================================
# Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
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


package_declare(damage_non_local
  DESCRIPTION "Package for Non-local damage constitutives laws Akantu")

package_declare_sources(damage_non_local
  model/solid_mechanics/materials/material_damage/material_damage_non_local.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.cc
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.hh
  model/solid_mechanics/materials/material_damage/material_mazars_non_local_tmpl.hh
  model/solid_mechanics/materials/material_damage/material_von_mises_mazars_non_local.cc
  model/solid_mechanics/materials/material_damage/material_von_mises_mazars_non_local.hh
  
  model/solid_mechanics/materials/weight_functions/damaged_weight_function.hh
  model/solid_mechanics/materials/weight_functions/damaged_weight_function.cc
  model/solid_mechanics/materials/weight_functions/damaged_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_weight_function.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_weight_function.cc
  model/solid_mechanics/materials/weight_functions/remove_damaged_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_with_damage_rate_weight_function.hh
  model/solid_mechanics/materials/weight_functions/remove_damaged_with_damage_rate_weight_function_inline_impl.hh
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function.hh
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function.cc
  model/solid_mechanics/materials/weight_functions/stress_based_weight_function_inline_impl.hh
  )
