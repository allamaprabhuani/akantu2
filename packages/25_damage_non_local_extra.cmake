#===============================================================================
# @file   damage_non_local.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Fri Jun 15 13:48:37 2012
#
# @brief  package description for non-local materials
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

option(AKANTU_DAMAGE_NON_LOCAL_EXTRA "Package for Non-local damage constitutives laws Akantu" OFF)

add_internal_package_dependencies(damage_non_local_extra extra_materials)
add_external_package_dependencies(damage_non_local_extra damage_non_local)

set(AKANTU_DAMAGE_NON_LOCAL_EXTRA_FILES
  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.cc

  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local.hh
  model/solid_mechanics/materials/material_damage/material_brittle_non_local.hh
  model/solid_mechanics/materials/material_damage/material_damage_iterative_non_local.hh

  model/solid_mechanics/materials/material_damage/material_vreepeerlings_non_local_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_brittle_non_local_inline_impl.cc
  model/solid_mechanics/materials/material_damage/material_damage_iterative_non_local_inline_impl.cc

  model/solid_mechanics/materials/material_non_local_extra_includes.hh
  )

set(AKANTU_DAMAGE_NON_LOCAL_EXTRA_TESTS
  )

set(AKANTU_DAMAGE_NON_LOCAL_EXTRA_DOCUMENTATION "
This package activates the non local damage feature of AKANTU
")
