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