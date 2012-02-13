set(AKANTU_CORE ON CACHE INTERNAL "core package for Akantu" FORCE)
set(CORE_FILES
  # source files
  common/aka_common.cc
  common/aka_error.cc
  common/aka_extern.cc
  common/aka_static_memory.cc
  common/aka_memory.cc
  common/aka_vector.cc
  common/aka_math.cc
  fem/shape_lagrange.cc
  fem/shape_linked.cc
  fem/shape_cohesive.cc
  fem/integrator_gauss.cc
  fem/mesh.cc
  fem/fem.cc
  fem/element_class.cc
  fem/cohesive_element.cc
  fem/fem_template.cc
  model/model.cc
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_boundary.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/material.cc
  model/solid_mechanics/materials/material_non_local.cc
  model/solid_mechanics/materials/material_elastic.cc
  model/solid_mechanics/materials/material_elastic_caughey.cc
  model/solid_mechanics/materials/material_viscoelastic.cc
  model/solid_mechanics/materials/material_damage.cc
  model/solid_mechanics/materials/material_marigo.cc
  model/solid_mechanics/materials/material_mazars.cc
  model/solid_mechanics/materials/material_neohookean.cc
  model/solid_mechanics/materials/material_marigo_non_local.cc
  model/solid_mechanics/materials/material_mazars_non_local.cc
  model/solid_mechanics/materials/material_damage_linear.cc
  mesh_utils/mesh_io.cc
  mesh_utils/mesh_pbc.cc
  mesh_utils/mesh_io/mesh_io_msh.cc
  mesh_utils/mesh_io/mesh_io_msh_struct.cc
  mesh_utils/mesh_io/mesh_io_diana.cc
  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_utils.cc
  solver/sparse_matrix.cc
  solver/solver.cc
  synchronizer/synchronizer_registry.cc
  synchronizer/synchronizer.cc
  synchronizer/distributed_synchronizer.cc
  synchronizer/pbc_synchronizer.cc
  synchronizer/data_accessor.cc
  synchronizer/static_communicator.cc
  synchronizer/grid_synchronizer.cc
  synchronizer/dof_synchronizer.cc

  #header files

  mesh_utils/mesh_io/mesh_io_msh.hh
  mesh_utils/mesh_io/mesh_io_msh_struct.hh
  mesh_utils/mesh_io/mesh_io_diana.hh
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_partition.hh
  mesh_utils/mesh_io.hh
  mesh_utils/mesh_partition/mesh_partition_scotch.hh
  solver/sparse_matrix.hh
  solver/solver.hh
  synchronizer/synchronizer.hh
  synchronizer/synchronizer_registry.hh
  synchronizer/static_communicator_dummy.hh
  synchronizer/static_communicator_inline_impl.hh
  synchronizer/distributed_synchronizer.hh
  synchronizer/pbc_synchronizer.hh
  synchronizer/static_communicator.hh
  synchronizer/dof_synchronizer.hh
  synchronizer/real_static_communicator.hh
  synchronizer/data_accessor.hh
  synchronizer/communication_buffer.hh
  synchronizer/grid_synchronizer.hh
  common/aka_grid.hh
  common/aka_grid_tmpl.hh
  common/aka_types.hh
  common/aka_static_memory.hh
  common/aka_static_memory_tmpl.hh
  common/aka_memory.hh
  common/aka_math.hh
  common/aka_math_tmpl.hh
  common/aka_csr.hh
  common/aka_error.hh
  common/aka_common.hh
  common/aka_vector.hh
  common/aka_vector_tmpl.hh
  common/aka_types_expression.hh
  common/aka_circular_vector.hh
  fem/mesh.hh
  fem/fem.hh
  fem/by_element_type.hh
  fem/shape_functions.hh
  fem/shape_lagrange.hh
  fem/shape_cohesive.hh
  fem/fem_template.hh
  fem/integrator_gauss.hh
  fem/integrator.hh
  fem/element_class.hh
  fem/cohesive_element.hh
  fem/shape_linked.hh
  model/model.hh
  model/parser.hh
  model/parser_tmpl.hh
  model/structural_mechanics/structural_mechanics_model.hh
  model/integration_scheme/integration_scheme_2nd_order.hh
  model/integration_scheme/generalized_trapezoidal.hh
  model/integration_scheme/newmark-beta.hh
  model/integration_scheme/integration_scheme_1st_order.hh
  model/solid_mechanics/materials/material_damage.hh
  model/solid_mechanics/materials/material_marigo.hh
  model/solid_mechanics/materials/material_marigo_non_local.hh
  model/solid_mechanics/materials/material_mazars_non_local.hh
  model/solid_mechanics/materials/material_elastic_caughey.hh
  model/solid_mechanics/materials/material_viscoelastic.hh
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_non_local.hh
  model/solid_mechanics/materials/material_mazars.hh
  model/solid_mechanics/materials/material_damage_linear.hh
  model/solid_mechanics/materials/material_neohookean.hh
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_tmpl.hh
  model/solid_mechanics/material.hh
  model/heat_transfer/heat_transfer_model.hh

  #inline implementation files
  mesh_utils/mesh_utils_inline_impl.cc
  solver/sparse_matrix_inline_impl.cc
  solver/solver_inline_impl.cc
  synchronizer/dof_synchronizer_inline_impl.cc
  synchronizer/communication_buffer_inline_impl.cc
  common/aka_memory_inline_impl.cc
  common/aka_static_memory_inline_impl.cc
  common/aka_circular_vector_inline_impl.cc
  fem/integrator_gauss_inline_impl.cc
  fem/element_classes/element_class_triangle_3_inline_impl.cc
  fem/element_classes/element_class_segment_2_inline_impl.cc
  fem/element_classes/element_class_quadrangle_4_inline_impl.cc
  fem/element_classes/element_class_quadrangle_8_inline_impl.cc
  fem/element_classes/element_class_bernoulli_beam_2_inline_impl.cc
  fem/element_classes/element_class_hexahedron_8_inline_impl.cc
  fem/element_classes/element_class_triangle_6_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_10_inline_impl.cc
  fem/element_classes/element_class_segment_3_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_4_inline_impl.cc
  fem/shape_functions_inline_impl.cc
  fem/mesh_inline_impl.cc
  fem/element_class_inline_impl.cc
  fem/by_element_type_tmpl.hh
  fem/fem_inline_impl.cc
  fem/shape_linked_inline_impl.cc
  fem/shape_lagrange_inline_impl.cc
  model/model_inline_impl.cc
  model/integration_scheme/generalized_trapezoidal_inline_impl.cc
  model/integration_scheme/newmark-beta_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_inline_impl.cc
  model/solid_mechanics/materials/material_mazars_inline_impl.cc
  model/solid_mechanics/materials/material_elastic_inline_impl.cc
  model/solid_mechanics/materials/material_non_local_inline_impl.cc
  model/solid_mechanics/materials/material_neohookean_inline_impl.cc
  model/solid_mechanics/materials/material_marigo_inline_impl.cc
  model/solid_mechanics/materials/material_damage_linear_inline_impl.cc
  model/solid_mechanics/material_inline_impl.cc
  model/parser_inline_impl.cc
  )

set(CORE_DEB_DEPEND
  libboost-dev
  )
