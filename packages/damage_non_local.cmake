option(AKANTU_DAMAGE_NON_LOCAL "Package for Non-local damage constitutives laws Akantu" OFF)
set(AKANTU_DAMAGE_NON_LOCAL_FILES
  model/solid_mechanics/materials/weight_function.cc
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.cc

  model/solid_mechanics/materials/weight_function.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.hh
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.hh

  model/solid_mechanics/materials/material_non_local.hh
  model/solid_mechanics/materials/material_non_local_inline_impl.cc
  )