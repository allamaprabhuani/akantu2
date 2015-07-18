#===============================================================================
# @file   package.cmake
#
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
#
# @brief  package description for interface-enriched generalized IGFEM
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
#===============================================================================

package_declare(IGFEM
  DESCRIPTION "Use Interface-enriched generalized FEM"
  DEPENDS CGAL)

package_declare_sources(igfem
  element_class_igfem.cc
  element_class_igfem.hh
  element_classes_igfem/element_class_igfem_segment_3_inline_impl.cc
  element_classes_igfem/element_class_igfem_triangle_4_inline_impl.cc
  element_classes_igfem/element_class_igfem_triangle_5_inline_impl.cc
  integrator_gauss_igfem.hh
  integrator_gauss_igfem_inline_impl.cc
  interpolation_element_igfem.cc
  interpolation_element_igfem_tmpl.hh
  geometrical_element_igfem.cc
  shape_igfem.hh
  shape_igfem.cc
  shape_igfem_inline_impl.cc
  fe_engine_template_tmpl_igfem.hh
  dumper_igfem_connectivity.hh
  dumper_igfem_elemental_field.hh
  dumper_igfem_generic_elemental_field.hh
  dumper_igfem_element_iterator.hh
  dumper_igfem_material_internal_field.hh
  dumper_igfem_quadrature_points_field.hh
  dumper_igfem_element_partition.hh

  igfem_helper.hh
  igfem_helper.cc
  igfem_enrichment.hh
  igfem_enrichment.cc
  igfem_enrichment_inline_impl.cc
  solid_mechanics_model_igfem.hh
  solid_mechanics_model_igfem.cc
  mesh_igfem_spherical_growing_gel.hh

  material_igfem/material_igfem_includes.hh
  material_igfem/material_igfem.hh
  material_igfem/material_igfem.cc
  material_igfem/material_igfem_inline_impl.cc
  material_igfem/material_igfem_elastic.hh
  material_igfem/material_igfem_elastic.cc
  material_igfem/material_igfem_elastic_inline_impl.cc

  material_igfem/igfem_internal_field.hh
  material_igfem/igfem_internal_field_tmpl.hh
  )

package_declare_material_infos(igfem
  LIST AKANTU_IGFEM_MATERIAL_LIST
  INCLUDE material_igfem_includes.hh
  )
