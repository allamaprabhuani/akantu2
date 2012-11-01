add_meta_package(PARALLEL "Add parallel support in Akantu" OFF MPI SCOTCH)

set(AKANTU_PARALLEL_TESTS
  test_solid_mechanics_model_bar_traction2d_parallel
  test_solid_mechanics_model_segment_parallel
  test_solid_mechanics_model_pbc_parallel
  test_synchronizer_communication
  test_dof_synchronizer
  )