#===============================================================================
# @file   CMakeLists.txt
# @author Richart Nicolas <nicolas.richart@epfl.ch>
# @date   Fri Sep 29 16:46:30 2010 
#
# @brief package description for cohesive elements
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
option(AKANTU_COHESIVE_ELEMENT "Use cohesive_element package of Akantu" OFF)
set(AKANTU_COHESIVE_ELEMENT_FILES
  model/solid_mechanics/materials/material_cohesive_includes.hh

  fem/shape_cohesive.cc
  fem/cohesive_element.cc
  fem/shape_cohesive.hh
  fem/cohesive_element.hh
  fem/fem_template_cohesive.cc

  fem/integrator_cohesive.hh
  fem/integrator_cohesive_inline_impl.cc
  fem/fem_template_inline_impl.cc
  fem/shape_cohesive_inline_impl.cc

  common/aka_common_inline_impl.cc
  model/solid_mechanics/materials/material_cohesive/material_cohesive_inline_impl.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_exponential_inline_impl.cc

  model/solid_mechanics/solid_mechanics_model_cohesive.cc
  model/solid_mechanics/materials/material_cohesive/material_cohesive.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_bilinear.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear_extrinsic.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_exponential.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear_exponential_extrinsic.cc

  model/solid_mechanics/solid_mechanics_model_cohesive.hh
  model/solid_mechanics/materials/material_cohesive/material_cohesive.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_bilinear.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_exponential.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear_extrinsic.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear_exponential_extrinsic.hh
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh
  )


set(AKANTU_COHESIVE_ELEMENT_TESTS
  )

set(AKANTU_COHESIVE_ELEMENT_DOC
  manual/manual-cohesive_element.tex
  )