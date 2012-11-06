/**
 * @file   material_cohesive_includes.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Oct 31 16:24:42 2012
 *
 * @brief  List of includes for cohesive elements
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

#include "material_cohesive.hh"
#include "material_cohesive_linear.hh"
#include "material_cohesive_bilinear.hh"
#include "material_cohesive_linear_extrinsic.hh"
#include "material_cohesive_exponential.hh"
#include "material_cohesive_linear_exponential_extrinsic.hh"

#define AKANTU_COHESIVE_MATERIAL_LIST					\
  ((2, (cohesive_bilinear      , MaterialCohesiveBilinear     )))	\
  ((2, (cohesive_linear        , MaterialCohesiveLinear       )))	\
  ((2, (cohesive_linear_extrinsic, MaterialCohesiveLinearExtrinsic )))	\
  ((2, (cohesive_linear_exponential_extrinsic, MaterialCohesiveLinearExponentialExtrinsic ))) \
  ((2, (cohesive_exponential   , MaterialCohesiveExponential  )))
