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
  common/aka_blas_lapack.hh
  common/aka_circular_vector.hh
  common/aka_circular_vector_inline_impl.cc
  common/aka_common.cc
  common/aka_common.hh
  common/aka_common_inline_impl.cc
  common/aka_csr.hh
  common/aka_error.cc
  common/aka_error.hh
  common/aka_event_handler.hh
  common/aka_extern.cc
  common/aka_fwd.hh
  common/aka_grid.hh
  common/aka_grid_dynamic.hh
  common/aka_grid_tmpl.hh
  common/aka_math.cc
  common/aka_math.hh
  common/aka_math_tmpl.hh
  common/aka_memory.cc
  common/aka_memory.hh
  common/aka_memory_inline_impl.cc
  common/aka_random_generator.hh
  common/aka_safe_enum.hh
  common/aka_static_memory.cc
  common/aka_static_memory.hh
  common/aka_static_memory_inline_impl.cc
  common/aka_static_memory_tmpl.hh
  common/aka_typelist.hh
  common/aka_types.hh
  common/aka_vector.cc
  common/aka_vector.hh
  common/aka_vector_tmpl.hh
  common/aka_visitor.hh
  common/aka_voigthelper.hh
  common/aka_voigthelper.cc

  fem/by_element_type.hh
  fem/by_element_type_tmpl.hh

  fem/element_class.cc
  fem/element_class.hh
  fem/element_class_tmpl.hh
  fem/element_classes/element_class_hexahedron_8_inline_impl.cc
  fem/element_classes/element_class_pentahedron_6_inline_impl.cc
  fem/element_classes/element_class_point_1_inline_impl.cc
  fem/element_classes/element_class_quadrangle_4_inline_impl.cc
  fem/element_classes/element_class_quadrangle_8_inline_impl.cc
  fem/element_classes/element_class_segment_2_inline_impl.cc
  fem/element_classes/element_class_segment_3_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_10_inline_impl.cc
  fem/element_classes/element_class_tetrahedron_4_inline_impl.cc
  fem/element_classes/element_class_triangle_3_inline_impl.cc
  fem/element_classes/element_class_triangle_6_inline_impl.cc

  fem/element_group.cc
  fem/element_group.hh
  fem/element_group_inline_impl.cc
  fem/fem.cc
  fem/fem.hh
  fem/fem_inline_impl.cc
  fem/fem_template.hh
  fem/fem_template_tmpl.hh
  fem/geometrical_data_tmpl.hh
  fem/geometrical_element.cc
  fem/group_manager.cc
  fem/group_manager.hh
  fem/group_manager_inline_impl.cc
  fem/integration_element.cc
  fem/integrator.hh
  fem/integrator_gauss.hh
  fem/integrator_gauss_inline_impl.cc
  fem/interpolation_element.cc
  fem/interpolation_element_tmpl.hh
  fem/mesh.cc
  fem/mesh.hh
  fem/mesh_filter.hh
  fem/mesh_data.cc
  fem/mesh_data.hh
  fem/mesh_data_tmpl.hh
  fem/mesh_inline_impl.cc
  fem/node_group.cc
  fem/node_group.hh
  fem/node_group_inline_impl.cc
  fem/shape_functions.hh
  fem/shape_functions_inline_impl.cc
  fem/shape_lagrange.cc
  fem/shape_lagrange.hh
  fem/shape_lagrange_inline_impl.cc
  fem/shape_linked.cc
  fem/shape_linked.hh
  fem/shape_linked_inline_impl.cc

  io/dumper/dumpable.hh
  io/dumper/dumpable_inline_impl.hh

  io/mesh_io.cc
  io/mesh_io.hh
  io/mesh_io/mesh_io_diana.cc
  io/mesh_io/mesh_io_diana.hh
  io/mesh_io/mesh_io_msh.cc
  io/mesh_io/mesh_io_msh.hh
  io/model_io.cc
  io/model_io.hh

  io/parser/algebraic_parser.hh
  io/parser/input_file_parser.hh
  io/parser/parsable.cc
  io/parser/parsable.hh
  io/parser/parsable_tmpl.hh
  io/parser/parser.cc
  io/parser/parser.hh
  io/parser/parser_tmpl.hh

  mesh_utils/mesh_partition.cc
  mesh_utils/mesh_partition.hh
  mesh_utils/mesh_partition/mesh_partition_mesh_data.cc
  mesh_utils/mesh_partition/mesh_partition_mesh_data.hh
  mesh_utils/mesh_partition/mesh_partition_scotch.hh
  mesh_utils/mesh_pbc.cc
  mesh_utils/mesh_utils.cc
  mesh_utils/mesh_utils.hh
  mesh_utils/mesh_utils_inline_impl.cc

  model/boundary_condition.hh
  model/boundary_condition_functor.hh
  model/boundary_condition_functor_inline_impl.cc
  model/boundary_condition_tmpl.hh
  model/integration_scheme/generalized_trapezoidal.hh
  model/integration_scheme/generalized_trapezoidal_inline_impl.cc
  model/integration_scheme/integration_scheme_1st_order.hh
  model/integration_scheme/integration_scheme_2nd_order.hh
  model/integration_scheme/newmark-beta.hh
  model/integration_scheme/newmark-beta_inline_impl.cc
  model/model.cc
  model/model.hh
  model/model_inline_impl.cc

  model/solid_mechanics/material.cc
  model/solid_mechanics/material.hh
  model/solid_mechanics/material_inline_impl.cc
  model/solid_mechanics/material_list.hh
  model/solid_mechanics/material_random_internal.hh
  model/solid_mechanics/material_selector.hh
  model/solid_mechanics/material_selector_tmpl.hh
  model/solid_mechanics/materials/internal_field.hh
  model/solid_mechanics/materials/internal_field_tmpl.hh
  model/solid_mechanics/materials/material_elastic.cc
  model/solid_mechanics/materials/material_elastic.hh
  model/solid_mechanics/materials/material_elastic_inline_impl.cc
  model/solid_mechanics/materials/material_thermal.cc
  model/solid_mechanics/materials/material_thermal.hh
  model/solid_mechanics/materials/random_internal_field.hh
  model/solid_mechanics/materials/random_internal_field_tmpl.hh
  model/solid_mechanics/solid_mechanics_model.cc
  model/solid_mechanics/solid_mechanics_model.hh
  model/solid_mechanics/solid_mechanics_model_inline_impl.cc
  model/solid_mechanics/solid_mechanics_model_mass.cc
  model/solid_mechanics/solid_mechanics_model_material.cc
  model/solid_mechanics/solid_mechanics_model_tmpl.hh

  solver/solver.cc
  solver/solver.hh
  solver/sparse_matrix.cc
  solver/sparse_matrix.hh
  solver/sparse_matrix_inline_impl.cc

  synchronizer/communication_buffer.hh
  synchronizer/communication_buffer_inline_impl.cc
  synchronizer/data_accessor.cc
  synchronizer/data_accessor.hh
  synchronizer/data_accessor_inline_impl.cc
  synchronizer/distributed_synchronizer.cc
  synchronizer/distributed_synchronizer.hh
  synchronizer/distributed_synchronizer_tmpl.hh
  synchronizer/dof_synchronizer.cc
  synchronizer/dof_synchronizer.hh
  synchronizer/dof_synchronizer_inline_impl.cc
  synchronizer/filtered_synchronizer.cc
  synchronizer/filtered_synchronizer.hh
  synchronizer/mpi_type_wrapper.hh
  synchronizer/pbc_synchronizer.cc
  synchronizer/pbc_synchronizer.hh
  synchronizer/real_static_communicator.hh
  synchronizer/static_communicator.cc
  synchronizer/static_communicator.hh
  synchronizer/static_communicator_dummy.hh
  synchronizer/static_communicator_inline_impl.hh
  synchronizer/synchronizer.cc
  synchronizer/synchronizer.hh
  synchronizer/synchronizer_registry.cc
  synchronizer/synchronizer_registry.hh
  )


set(AKANTU_CORE_DEB_DEPEND
  libboost-dev
  )

set(AKANTU_CORE_TESTS
  test_csr
  test_facet_element_mapping
  test_facet_extraction_tetrahedron_4
  test_facet_extraction_triangle_3
  test_grid
  test_interpolate_stress
  test_local_material
  test_material_damage_non_local
  test_material_thermal
  test_matrix
  test_mesh_boundary
  test_mesh_data
  test_mesh_io_msh
  test_mesh_io_msh_physical_names
  test_mesh_partitionate_mesh_data
  test_parser
  test_pbc_tweak
  test_purify_mesh
  test_solid_mechanics_model_bar_traction2d
  test_solid_mechanics_model_bar_traction2d_structured
  test_solid_mechanics_model_bar_traction2d_structured_pbc
  test_solid_mechanics_model_boundary_condition
  test_solid_mechanics_model_circle_2
  test_solid_mechanics_model_cube3d
  test_solid_mechanics_model_cube3d_pbc
  test_solid_mechanics_model_cube3d_tetra10
  test_solid_mechanics_model_square
  test_solid_mechanics_model_reassign_material
  test_solid_mechanics_model_prestrain_material
  test_static_memory
  test_surface_extraction_tetrahedron_4
  test_surface_extraction_triangle_3
  test_vector
  test_vector_iterator
  test_weight
  )

set(AKANTU_CORE_MANUAL_FILES
  manual.sty
  manual.cls
  manual.tex
  manual-macros.sty
  manual-titlepages.tex
  manual-introduction.tex
  manual-gettingstarted.tex
  manual-io.tex
  manual-solidmechanicsmodel.tex
  manual-lumping.tex
  manual-elements.tex
  manual-appendix-elements.tex
  manual-backmatter.tex
  manual-bibliography.bib
  manual-bibliographystyle.bst

  figures/bc_and_ic_example.pdf
  figures/boundary.pdf
  figures/boundary.svg
  figures/dirichlet.pdf
  figures/dirichlet.svg
  figures/doc_wheel.pdf
  figures/doc_wheel.svg
  figures/dynamic_analysis.png
  figures/explicit_dynamic.pdf
  figures/explicit_dynamic.svg
  figures/hooke_law.pdf
  figures/hot-point-1.png
  figures/hot-point-2.png
  figures/implicit_dynamic.pdf
  figures/implicit_dynamic.svg
  figures/implicit_static.pdf
  figures/implicit_static.svg
  figures/insertion.pdf
  figures/interpolate.pdf
  figures/interpolate.svg
  figures/law.pdf
  figures/static_analysis.png
  figures/stress_strain_el.pdf
  figures/tangent.pdf
  figures/tangent.svg
  figures/vectors.pdf
  figures/vectors.svg


  figures/elements/hexahedron_8.pdf
  figures/elements/hexahedron_8.svg
  figures/elements/quadrangle_4.pdf
  figures/elements/quadrangle_4.svg
  figures/elements/quadrangle_8.pdf
  figures/elements/quadrangle_8.svg
  figures/elements/segment_2.pdf
  figures/elements/segment_2.svg
  figures/elements/segment_3.pdf
  figures/elements/segment_3.svg
  figures/elements/tetrahedron_10.pdf
  figures/elements/tetrahedron_10.svg
  figures/elements/tetrahedron_4.pdf
  figures/elements/tetrahedron_4.svg
  figures/elements/triangle_3.pdf
  figures/elements/triangle_3.svg
  figures/elements/triangle_6.pdf
  figures/elements/triangle_6.svg
  figures/elements/xtemp.pdf
  )

find_program(READLINK_COMMAND readlink)
find_program(ADDR2LINE_COMMAND addr2line)
mark_as_advanced(READLINK_COMMAND)
mark_as_advanced(ADDR2LINE_COMMAND)

include(CheckFunctionExists)

check_function_exists(clock_gettime _clock_gettime)

if(NOT _clock_gettime)
  set(AKANTU_USE_OBSOLETE_GETTIMEOFDAY ON)
else()
  set(AKANTU_USE_OBSOLETE_GETTIMEOFDAY OFF)
endif()

