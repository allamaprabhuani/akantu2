#===============================================================================
# @file   core.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for core
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

set(AKANTU_CORE ON CACHE INTERNAL "core package for Akantu" FORCE)

set(AKANTU_CORE_FILES
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
  fem/mesh.cc
  fem/fem.cc
  fem/mesh_data.cc
  fem/boundary.cc
  fem/sub_boundary.cc
  io/dumper/dumpable.hh
  io/mesh_io.cc
  io/model_io.cc
  io/mesh_io/mesh_io_msh.cc
  io/mesh_io/mesh_io_diana.cc
  model/model.cc
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_boundary.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/material.cc
  model/solid_mechanics/material_parameters.cc
  model/solid_mechanics/materials/material_elastic.cc
  mesh_utils/mesh_pbc.cc
  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_partition/mesh_partition_mesh_data.cc
  mesh_utils/mesh_utils.cc
  solver/sparse_matrix.cc
  solver/solver.cc
  synchronizer/synchronizer_registry.cc
  synchronizer/synchronizer.cc
  synchronizer/distributed_synchronizer.cc
  synchronizer/pbc_synchronizer.cc
  synchronizer/data_accessor.cc
  synchronizer/static_communicator.cc
  synchronizer/dof_synchronizer.cc

  #header files
  io/mesh_io.hh
  io/model_io.hh
  io/mesh_io/mesh_io_msh.hh
  io/mesh_io/mesh_io_diana.hh
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_partition.hh
  mesh_utils/mesh_partition/mesh_partition_scotch.hh
  mesh_utils/mesh_partition/mesh_partition_mesh_data.hh
  solver/sparse_matrix.hh
  solver/solver.hh
  synchronizer/synchronizer.hh
  synchronizer/synchronizer_registry.hh
  synchronizer/static_communicator_dummy.hh
  synchronizer/static_communicator_inline_impl.hh
  synchronizer/distributed_synchronizer.hh
  synchronizer/distributed_synchronizer_tmpl.hh
  synchronizer/pbc_synchronizer.hh
  synchronizer/static_communicator.hh
  synchronizer/dof_synchronizer.hh
  synchronizer/real_static_communicator.hh
  synchronizer/data_accessor.hh
  synchronizer/communication_buffer.hh
  common/aka_fwd.hh
  common/aka_grid.hh
  common/aka_grid_tmpl.hh
  common/aka_types.hh
  common/aka_static_memory.hh
  common/aka_static_memory_tmpl.hh
  common/aka_memory.hh
  common/aka_math.hh
  common/aka_math_tmpl.hh
  common/aka_blas_lapack.hh
  common/aka_csr.hh
  common/aka_error.hh
  common/aka_common.hh
  common/aka_vector.hh
  common/aka_vector_tmpl.hh
  common/aka_circular_vector.hh
  common/aka_event_handler.hh
  common/aka_random_generator.hh
  common/aka_bounding_box.hh
  common/aka_ci_string.hh
  common/aka_plane.hh
  common/aka_polytope.hh
  common/aka_sphere.hh
  common/aka_timer.hh
  common/aka_tree.hh
  common/aka_typelist.hh
  common/aka_visitor.hh
  common/aka_grid_dynamic.hh
  common/aka_safe_enum.hh
  fem/mesh.hh
  fem/fem.hh
  fem/by_element_type.hh
  fem/shape_functions.hh
  fem/shape_lagrange.hh
  fem/fem_template.hh
  fem/fem_template_tmpl.hh
  fem/integrator_gauss.hh
  fem/integrator.hh
  fem/element_class.hh
  fem/shape_linked.hh
  fem/geometrical_data_tmpl.hh
  fem/mesh_data.hh
  fem/boundary.hh
  fem/sub_boundary.hh
  model/model.hh
  model/boundary_condition.hh
  model/boundary_condition_functor.hh
  model/parser.hh
  model/parser_tmpl.hh
  model/structural_mechanics/structural_mechanics_model.hh
  model/integration_scheme/integration_scheme_2nd_order.hh
  model/integration_scheme/generalized_trapezoidal.hh
  model/integration_scheme/newmark-beta.hh
  model/integration_scheme/integration_scheme_1st_order.hh
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_tmpl.hh
  model/solid_mechanics/material.hh
  model/solid_mechanics/material_parameters.hh
  model/solid_mechanics/material_parameters_tmpl.hh
  model/solid_mechanics/materials/material_elastic.hh

  #inline implementation files
  mesh_utils/mesh_utils_inline_impl.cc
  solver/sparse_matrix_inline_impl.cc
  synchronizer/dof_synchronizer_inline_impl.cc
  synchronizer/communication_buffer_inline_impl.cc
  common/aka_common_inline_impl.cc
  common/aka_memory_inline_impl.cc
  common/aka_static_memory_inline_impl.cc
  common/aka_circular_vector_inline_impl.cc
  fem/integrator_gauss_inline_impl.cc
  fem/element_classes/element_class_point_1_inline_impl.cc
  fem/element_classes/element_class_triangle_3_inline_impl.cc
  fem/element_classes/element_class_segment_2_inline_impl.cc
  fem/element_classes/element_class_quadrangle_4_inline_impl.cc
  fem/element_classes/element_class_quadrangle_8_inline_impl.cc
  fem/element_classes/element_class_hexahedron_8_inline_impl.cc
  fem/element_classes/element_class_triangle_6_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_10_inline_impl.cc
  fem/element_classes/element_class_segment_3_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_4_inline_impl.cc
  fem/shape_functions_inline_impl.cc
  fem/mesh_inline_impl.cc
  fem/boundary_inline_impl.cc
  fem/sub_boundary_inline_impl.cc
  fem/mesh_data_tmpl.hh
  fem/by_element_type_tmpl.hh
  fem/fem_inline_impl.cc
  fem/shape_linked_inline_impl.cc
  fem/shape_lagrange_inline_impl.cc
  model/model_inline_impl.cc
  model/boundary_condition_functor_inline_impl.cc
  model/integration_scheme/generalized_trapezoidal_inline_impl.cc
  model/integration_scheme/newmark-beta_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_inline_impl.cc
  model/solid_mechanics/materials/material_elastic_inline_impl.cc
  model/solid_mechanics/material_inline_impl.cc
  model/parser_inline_impl.cc
  model/boundary_condition_tmpl.hh
  fem/geometrical_element.cc
  fem/element_class_tmpl.hh
  fem/element_class.cc
  fem/integration_element.cc
  fem/interpolation_element.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.cc
  model/solid_mechanics/materials/material_elastic_orthotropic.hh
  model/solid_mechanics/materials/material_elastic_orthotropic_inline_impl.cc
  )


set(AKANTU_CORE_DEB_DEPEND
  libboost-dev
  )

set(AKANTU_CORE_TESTS
  test_solid_mechanics_model_square
  test_vector
  test_vector_iterator
  test_matrix
  test_csr
  test_grid
  test_static_memory
  test_mesh_io_msh
  test_mesh_io_msh_physical_names
  test_mesh_data
  test_facet_element_mapping
  test_mesh_boundary
  test_facet_extraction_triangle_3
  test_facet_extraction_tetrahedron_4
  test_pbc_tweak
  test_purify_mesh
  test_local_material
  test_interpolate_stress
  test_weight
  test_solid_mechanics_model_boundary_condition
  test_solid_mechanics_model_circle_2
  test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d_structured
  test_solid_mechanics_model_bar_traction2d_structured_pbc
  test_solid_mechanics_model_cube3d
  test_solid_mechanics_model_cube3d_tetra10
  test_solid_mechanics_model_cube3d_pbc
  test_surface_extraction_triangle_3
  test_surface_extraction_tetrahedron_4
  test_material_damage_non_local
  )

find_program(READLINK_COMMAND readlink)
find_program(ADDR2LINE_COMMAND addr2line)
mark_as_advanced(READLINK_COMMAND)
mark_as_advanced(ADDR2LINE_COMMAND)

