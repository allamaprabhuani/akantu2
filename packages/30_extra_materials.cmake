#===============================================================================
# @file   extra_materials.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Wed Oct 31 16:24:42 2012
#
# @brief  package description for extra materials list
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
#===============================================================================

option(AKANTU_EXTRA_MATERIALS "Add the extra list of materials in Akantu" OFF)
add_external_package_dependencies(extra_materials lapack)

set(AKANTU_EXTRA_MATERIALS_FILES
  model/solid_mechanics/materials/material_extra_includes.hh

  model/solid_mechanics/materials/material_damage/material_brittle.cc
  model/solid_mechanics/materials/material_damage/material_brittle.hh
  model/solid_mechanics/materials/material_damage/material_brittle_inline_impl.cc

  model/solid_mechanics/materials/material_damage/material_damage.hh
  model/solid_mechanics/materials/material_damage/material_damage_tmpl.hh

  model/solid_mechanics/materials/material_damage/material_damage_iterative.cc
  model/solid_mechanics/materials/material_damage/material_damage_iterative.hh
  model/solid_mechanics/materials/material_damage/material_damage_iterative_inline_impl.cc

  model/solid_mechanics/materials/material_damage/material_damage_linear.cc
  model/solid_mechanics/materials/material_damage/material_damage_linear.hh
  model/solid_mechanics/materials/material_damage/material_damage_linear_inline_impl.cc

  model/solid_mechanics/materials/material_damage/material_marigo.cc
  model/solid_mechanics/materials/material_damage/material_marigo.hh
  model/solid_mechanics/materials/material_damage/material_marigo_inline_impl.cc

  model/solid_mechanics/materials/material_damage/material_mazars.cc
  model/solid_mechanics/materials/material_damage/material_mazars.hh
  model/solid_mechanics/materials/material_damage/material_mazars_inline_impl.cc

  model/solid_mechanics/materials/material_damage/material_vreepeerlings.hh
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_tmpl.hh

  model/solid_mechanics/materials/material_elastic_linear_anisotropic.cc
  model/solid_mechanics/materials/material_elastic_linear_anisotropic.hh

  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.hh

  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.cc
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean.hh
  model/solid_mechanics/materials/material_finite_deformation/material_neohookean_inline_impl.cc

  model/solid_mechanics/materials/material_plastic/material_plasticityinc.cc
  model/solid_mechanics/materials/material_plastic/material_plasticityinc.hh
  model/solid_mechanics/materials/material_plastic/material_plasticityinc_inline_impl.cc

  model/solid_mechanics/materials/material_plastic/material_viscoplasticity.cc
  model/solid_mechanics/materials/material_plastic/material_viscoplasticity.hh
  model/solid_mechanics/materials/material_plastic/material_viscoplasticity_inline_impl.cc

  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh

  model/solid_mechanics/materials/material_viscoelastic/material_stiffness_proportional.cc
  model/solid_mechanics/materials/material_viscoelastic/material_stiffness_proportional.hh

  )

set(AKANTU_EXTRA_MATERIALS_TESTS
  test_material_standard_linear_solid_deviatoric_relaxation
  test_material_standard_linear_solid_deviatoric_relaxation_tension
  test_material_plasticity
  )

set(AKANTU_EXTRA_MATERIALS_DOC
  manual/manual-extra_materials.tex
  )

set(AKANTU_EXTRA_MATERIALS_MANUAL_FILES
  manual-extra_materials.tex

  figures/stress_strain_neo.pdf
  figures/stress_strain_visco.pdf
  figures/visco_elastic_law.pdf
  )