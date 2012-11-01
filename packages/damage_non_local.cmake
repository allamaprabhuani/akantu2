option(AKANTU_DAMAGE_NON_LOCAL "Package for Non-local damage constitutives laws Akantu" OFF)

add_package_dependecies(damage_non_local extra_materials)

set(AKANTU_DAMAGE_NON_LOCAL_FILES
  model/solid_mechanics/materials/weight_function.cc
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.cc
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.cc

  model/solid_mechanics/materials/material_non_local_includes.hh
  model/solid_mechanics/materials/weight_function.hh
  model/solid_mechanics/materials/material_damage/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_damage/material_mazars_non_local.hh
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.hh

  model/solid_mechanics/materials/material_non_local.hh
  model/solid_mechanics/materials/material_non_local_inline_impl.cc
  model/solid_mechanics/materials/weight_function_tmpl.hh

  model/solid_mechanics/materials/material_damage/material_marigo_non_local_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_damage_non_local.hh
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local_inline_impl.cc

  synchronizer/grid_synchronizer.cc
  synchronizer/grid_synchronizer.hh
  )

set(AKANTU_DAMAGE_NON_LOCAL_TESTS
  test_material_damage_non_local
  test_grid_synchronizer
  )