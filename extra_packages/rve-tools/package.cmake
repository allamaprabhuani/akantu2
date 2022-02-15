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

package_declare(rve_tools
  DESCRIPTION "Multi-scale model package"
  DEPENDS extra_materials cohesive_element)

package_declare_sources(rve_tools
  rve_tools.cc
  rve_tools.hh
  material_FE2/material_FE2.hh
  material_FE2/material_FE2.cc
  solid_mechanics_model_RVE.hh
  solid_mechanics_model_RVE.cc
  mesh_utils/nodes_flag_updater.hh
  mesh_utils/nodes_flag_updater.cc
  mesh_utils/nodes_flag_updater_inline_impl.hh
  mesh_utils/nodes_eff_stress_updater.hh
  mesh_utils/nodes_eff_stress_updater.cc
  mesh_utils/nodes_eff_stress_updater_inline_impl.hh
  mesh_utils/crack_numbers_updater.hh
  mesh_utils/crack_numbers_updater.cc
  mesh_utils/crack_numbers_updater_inline_impl.hh
)
