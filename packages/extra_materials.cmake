option(AKANTU_EXTRA_MATERIALS "Add the extra list of materials in Akantu" OFF)

set(AKANTU_EXTRA_MATERIALS_FILES
  model/solid_mechanics/materials/material_neohookean.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_viscoelastic/material_stiffness_proportional.cc
  model/solid_mechanics/materials/material_damage/material_damage.cc
  model/solid_mechanics/materials/material_damage/material_marigo.cc
  model/solid_mechanics/materials/material_damage/material_mazars.cc
  model/solid_mechanics/materials/material_damage/material_damage_linear.cc
  model/solid_mechanics/materials/material_damage/material_vreepeerlings.cc

  model/solid_mechanics/materials/material_extra_includes.hh
  model/solid_mechanics/materials/material_viscoelastic/material_stiffness_proportional.hh
  model/solid_mechanics/materials/material_elastic_orthotropic.hh
  model/solid_mechanics/materials/material_neohookean.hh
  model/solid_mechanics/materials/material_damage/material_damage.hh
  model/solid_mechanics/materials/material_damage/material_marigo.hh
  model/solid_mechanics/materials/material_damage/material_mazars.hh
  model/solid_mechanics/materials/material_damage/material_damage_linear.hh
  model/solid_mechanics/materials/material_damage/material_vreepeerlings.hh

  model/solid_mechanics/materials/material_elastic_orthotropic_inline_impl.cc
  model/solid_mechanics/materials/material_neohookean_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_damage_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_marigo_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_mazars_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_damage_linear_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_inline_impl.cc

  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.cc
  model/solid_mechanics/materials/material_viscoelastic/material_standard_linear_solid_deviatoric.hh
  )

set(AKANTU_EXTRA_MATERIALS_TESTS
  test_material_standard_linear_solid_deviatoric_relaxation
  test_material_standard_linear_solid_deviatoric_relaxation_tension
  )

set(AKANTU_EXTRA_MATERIALS_DOC
  manual/manual-extra_materials.tex
  )