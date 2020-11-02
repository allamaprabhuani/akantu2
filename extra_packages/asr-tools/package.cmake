#===============================================================================
# @file   package.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @brief  package description for asr stuff
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
#===============================================================================

package_declare(asr_tools
  DESCRIPTION "ASR stuffs materials, model FE2, toolboxes"
  DEPENDS extra_materials)

package_declare_sources(asr_tools
  asr_tools.cc
  asr_tools.hh
  material_FE2/material_FE2.hh
  material_FE2/material_FE2.cc
  material_FE2/material_FE2_inline_impl.cc
  solid_mechanics_model_RVE.hh
  solid_mechanics_model_RVE.cc
  mesh_utils/nodes_pos_updater.hh
  mesh_utils/nodes_pos_updater.cc
  mesh_utils/nodes_pos_updater_inline_impl.cc
  mesh_utils/nodes_flag_updater.hh
  mesh_utils/nodes_flag_updater.cc
  mesh_utils/nodes_flag_updater_inline_impl.cc
  mesh_utils/crack_numbers_updater.hh
  mesh_utils/crack_numbers_updater.cc
  mesh_utils/crack_numbers_updater_inline_impl.cc
)
