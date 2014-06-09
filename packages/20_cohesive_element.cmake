#===============================================================================
# @file   cohesive_element.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  package description for cohesive elements
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

option(AKANTU_COHESIVE_ELEMENT "Use cohesive_element package of Akantu" OFF)
add_external_package_dependencies(cohesive_element lapack)

set(AKANTU_COHESIVE_ELEMENT_FILES
  model/solid_mechanics/materials/material_cohesive_includes.hh

  mesh_utils/cohesive_element_inserter.hh
  mesh_utils/cohesive_element_inserter.cc

  fem/cohesive_element.cc
  fem/shape_cohesive.hh
  fem/cohesive_element.hh
  fem/fe_engine_template_cohesive.cc

  fem/shape_cohesive_inline_impl.cc

  model/solid_mechanics/materials/material_cohesive/cohesive_internal_field_tmpl.hh
  model/solid_mechanics/materials/material_cohesive/cohesive_internal_field.hh

  model/solid_mechanics/materials/material_cohesive/material_cohesive_inline_impl.cc

  model/solid_mechanics/solid_mechanics_model_cohesive.cc
  model/solid_mechanics/fragment_manager.cc
  model/solid_mechanics/materials/material_cohesive/material_cohesive.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_bilinear.cc
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_exponential.cc

  model/solid_mechanics/solid_mechanics_model_cohesive.hh
  model/solid_mechanics/fragment_manager.hh
  model/solid_mechanics/materials/material_cohesive/material_cohesive.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_bilinear.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_linear.hh
  model/solid_mechanics/materials/material_cohesive/constitutive_laws/material_cohesive_exponential.hh

  io/dumper/dumper_iohelper_tmpl_connectivity_field.hh
  )


set(AKANTU_COHESIVE_ELEMENT_TESTS
  test_cohesive_buildfacets_tetrahedron
  test_cohesive_buildfacets_hexahedron
  test_cohesive_intrinsic
  test_cohesive_intrinsic_quadrangle
  test_cohesive_extrinsic
  test_cohesive_extrinsic_quadrangle
  test_cohesive_buildfragments
  test_cohesive_intrinsic_impl
  )

set(AKANTU_COHESIVE_ELEMENT_DOC
  manual/manual-cohesive_element.tex
  )

set(AKANTU_COHESIVE_ELEMENT_MANUAL_FILES
  manual-cohesive_element.tex

  figures/cohesive2d.pdf
  figures/cohesive3d.pdf
  figures/cohesive_exponential.pdf
  figures/linear_cohesive_law.pdf
  figures/bilinear_cohesive_law.pdf
  )
