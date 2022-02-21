#===============================================================================
# @file   package.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Mar 15 2018
# @date last modification:  Wed Dec 09 2020
#
# @brief  package description for extra materials list
#
#
# @section LICENSE
#
# Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(extra_materials ADVANCED
  DESCRIPTION "Add the extra list of materials in Akantu"
  DEPENDS lapack)

package_declare_sources(extra_materials
  material_extra_includes.hh

  material_damage/material_brittle.cc
  material_damage/material_brittle.hh
  material_damage/material_brittle_inline_impl.hh

  material_damage/material_damage_iterative.cc
  material_damage/material_damage_iterative.hh
  material_damage/material_damage_iterative_inline_impl.hh

  material_damage/material_damage_iterative_isotropic.cc
  material_damage/material_damage_iterative_isotropic.hh
  material_damage/material_damage_iterative_isotropic_inline_impl.hh

  material_damage/material_damage_iterative_orthotropic.cc
  material_damage/material_damage_iterative_orthotropic.hh
  material_damage/material_damage_iterative_orthotropic_inline_impl.hh

  material_damage/material_iterative_stiffness_reduction.cc
  material_damage/material_iterative_stiffness_reduction.hh
  material_damage/material_iterative_stiffness_reduction_inline_impl.hh

  material_damage/material_iterative_stiffness_reduction_orthotropic.cc
  material_damage/material_iterative_stiffness_reduction_orthotropic.hh

  material_damage/material_damage_iterative_viscoelastic.cc
  material_damage/material_damage_iterative_viscoelastic.hh
  material_damage/material_damage_iterative_viscoelastic_inline_impl.hh

  material_damage/material_damage_explicit.cc
  material_damage/material_damage_explicit.hh
  material_damage/material_damage_explicit_inline_impl.hh

  material_damage/material_damage_viscoelastic.cc
  material_damage/material_damage_viscoelastic.hh

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

  material_cohesive/material_cohesive_linear_sequential.hh
  material_cohesive/material_cohesive_linear_sequential.cc
  material_cohesive/material_cohesive_linear_sequential_inline_impl.hh
)

package_declare_material_infos(extra_materials
  LIST AKANTU_EXTRA_MATERIAL_LIST
  INCLUDE material_extra_includes.hh
  )

package_declare_documentation_files(extra_materials
  manual-extra_materials.tex
  manual-appendix-materials-extra-materials.tex

  figures/stress_strain_visco.pdf
  )

package_declare_documentation(extra_materials
  "This package activates additional constitutive laws:"
  "\\begin{itemize}"
  "  \\item Linear anisotropy"
  "  \\item Linear orthotropy"
  "  \\item Visco-plastic"
  "\\end{itemize}"
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

package_declare_material_infos(extra_materials_non_local
  LIST AKANTU_DAMAGE_NON_LOCAL_MATERIAL_EXTRA_LIST
  INCLUDE material_extra_includes.hh
  )
