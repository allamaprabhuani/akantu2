add_meta_package(IMPLICIT "Add support for implicit time scheme" OFF MUMPS MPI SCOTCH)

set(AKANTU_IMPLICIT_TESTS
  test_solid_mechanics_model_bar_traction2d_mass_not_lumped
  test_solid_mechanics_model_implicit_1d
  test_solid_mechanics_model_implicit_2d
  test_solid_mechanics_model_implicit_dynamic_2d
  )