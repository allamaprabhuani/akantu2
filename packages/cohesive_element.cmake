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


package_declare(cohesive_element
  DESCRIPTION "Use cohesive_element package of Akantu"
  DEPENDS lapack solid_mechanics)

package_declare_sources(cohesive_element
  fe_engine/cohesive_element.hh
  fe_engine/fe_engine_template_cohesive.cc
  fe_engine/shape_cohesive.hh
  fe_engine/shape_cohesive_inline_impl.hh

  mesh_utils/cohesive_element_inserter.cc
  mesh_utils/cohesive_element_inserter.hh
  mesh_utils/cohesive_element_inserter_inline_impl.hh
  mesh_utils/cohesive_element_inserter_parallel.cc
  mesh_utils/cohesive_element_inserter_helper.cc
  mesh_utils/cohesive_element_inserter_helper.hh

  model/solid_mechanics/solid_mechanics_model_cohesive/fragment_manager.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/fragment_manager.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/material_selector_cohesive.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/material_selector_cohesive.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/cohesive_internal_field.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/cohesive_internal_field_tmpl.hh

  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_bilinear.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_bilinear.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_exponential.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_exponential.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_fatigue.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_fatigue.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_friction.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_friction.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_uncoupled.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_uncoupled.hh

  model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive_includes.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive_inline_impl.hh

  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.cc
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive_inline_impl.hh
  model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive_parallel.cc
  )

package_declare_material_infos(cohesive_element
  LIST AKANTU_COHESIVE_MATERIAL_LIST
  INCLUDE material_cohesive_includes.hh
  )
