/**
 * @file   material_list.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Oct 19 2014
 * @date last modification:  Mon Feb 08 2021
 *
 * @brief  List of materials and all includes
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_LIST_HH__
#define __AKANTU_MATERIAL_LIST_HH__

#include "aka_config.hh"

/* -------------------------------------------------------------------------- */
/* Material includes                                                          */
/* -------------------------------------------------------------------------- */
#include "material_cohesive_includes.hh"
#include "material_core_includes.hh"
#include "material_non_local_includes.hh"

/* -------------------------------------------------------------------------- */
/* Material list                                                              */
/* -------------------------------------------------------------------------- */
#define AKANTU_MATERIAL_LIST                                                   \
  AKANTU_COHESIVE_MATERIAL_LIST                                                \
  AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST                                        \
  AKANTU_CORE_MATERIAL_LIST

// leave an empty line after the list

#endif /* __AKANTU_MATERIAL_LIST_HH__ */
