#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(extra_materials ADVANCED
  DESCRIPTION "Add the extra list of materials in Akantu"
  DEPENDS solid_mechanics)

package_declare_sources(extra_materials
  material_extra_includes.hh

  material_damage/material_brittle.cc
  material_damage/material_brittle.hh
  material_damage/material_brittle_inline_impl.hh

  material_damage/material_damage_iterative.cc
  material_damage/material_damage_iterative.hh
  material_damage/material_damage_iterative_inline_impl.hh

  material_damage/material_iterative_stiffness_reduction.cc
  material_damage/material_iterative_stiffness_reduction.hh

  material_damage/material_damage_linear.cc
  material_damage/material_damage_linear.hh
  material_damage/material_damage_linear_inline_impl.hh

  material_damage/material_vreepeerlings.hh
  material_damage/material_vreepeerlings_inline_impl.hh
  material_damage/material_vreepeerlings_tmpl.hh

  material_plastic/material_viscoplastic.cc
  material_plastic/material_viscoplastic.hh
  material_plastic/material_viscoplastic_inline_impl.hh

  material_viscoelastic/material_stiffness_proportional.cc
  material_viscoelastic/material_stiffness_proportional.hh

  material_damage/material_orthotropic_damage.hh
  material_damage/material_orthotropic_damage_tmpl.hh

  material_damage/material_orthotropic_damage_iterative.cc
  material_damage/material_orthotropic_damage_iterative.hh
  material_damage/material_orthotropic_damage_iterative_inline_impl.hh

  material_FE2/material_FE2.hh
  material_FE2/material_FE2.cc
  material_FE2/material_FE2_inline_impl.hh
  material_FE2/solid_mechanics_model_RVE.hh
  material_FE2/solid_mechanics_model_RVE.cc
  )

package_declare(extra_materials_non_local ADVANCED
  DESCRIPTION "Add the extra list of non local materials in Akantu"
  DEPENDS extra_materials damage_non_local)

package_declare_sources(extra_materials_non_local
  material_damage/material_orthotropic_damage_non_local.hh

  material_damage/material_vreepeerlings_non_local.cc
  material_damage/material_vreepeerlings_non_local.hh
  material_damage/material_brittle_non_local.hh
  material_damage/material_damage_iterative_non_local.hh
  material_damage/material_damage_iterative_non_local.cc
  material_damage/material_orthotropic_damage_iterative_non_local.hh

  material_damage/material_vreepeerlings_non_local_inline_impl.hh
  material_damage/material_brittle_non_local_inline_impl.hh
  material_damage/material_damage_iterative_non_local_inline_impl.hh
  material_damage/material_orthotropic_damage_iterative_non_local_inline_impl.hh

  material_non_local_extra_includes.hh
  )
