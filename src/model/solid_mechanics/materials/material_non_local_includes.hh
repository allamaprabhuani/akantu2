/**
 * @file   material_non_local_includes.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Oct 30 15:36:08 2012
 *
 * @brief  Non local materials includes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#  include "material_marigo_non_local.hh"
#  include "material_mazars_non_local.hh"
#  include "material_vreepeerlings_non_local.hh"

#define AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST			\
  ((stress_wf, StressBasedWeightFunction  ))				\
  ((damage_wf, DamagedWeightFunction      ))				\
  ((remove_wf, RemoveDamagedWeightFunction))				\
  ((base_wf,   BaseWeightFunction         ))


#define AKANTU_MATERIAL_VREEPEERLINGS_WEIGHT_FUNCTION_TMPL_LIST		\
  AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST				\
  ((removed_damrate_wf, RemoveDamagedWithDamageRateWeightFunction))

#define AKANTU_DAMAGE_NON_LOCAL_MATERIAL_LIST				\
  ((3, (marigo_non_local       , MaterialMarigoNonLocal,		\
	AKANTU_MATERIAL_WEIGHT_FUNCTION_TMPL_LIST)))			\
  ((2, (mazars_non_local       , MaterialMazarsNonLocal       )))	\
  ((3, (vreepeerlings_non_local, MaterialVreePeerlingsNonLocal,		\
	AKANTU_MATERIAL_VREEPEERLINGS_WEIGHT_FUNCTION_TMPL_LIST)))
