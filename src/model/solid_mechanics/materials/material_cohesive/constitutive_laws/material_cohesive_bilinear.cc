/**
 * @file   material_cohesive_bilinear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Feb 22 16:31:20 2012
 *
 * @brief  Bilinear cohesive constitutive law
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
#include "material_cohesive_bilinear.hh"
#include "solid_mechanics_model_cohesive.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::MaterialCohesiveBilinear(SolidMechanicsModel & model, const ID & id) :
  MaterialCohesiveLinear<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("delta_0", delta_0, 0.,
		      _pat_parsable | _pat_readable,
		      "Elastic limit displacement");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  /**
   * Recompute sigma_c as
   * @f$ {\sigma_c}_\textup{new} =
   * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
   */

  current_delta_c = 2 * this->G_cI / this->sigma_c;

  AKANTU_DEBUG_ASSERT(current_delta_c > delta_0, "Check your material file");

  Real sigma_c = current_delta_c / (current_delta_c - delta_0);

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  this->sigma_c_eff     .setDefaultValue(sigma_c);
  this->delta_c         .setDefaultValue(current_delta_c);
  this->delta_max       .setDefaultValue(delta_0);
  this->insertion_stress.setDefaultValue(0);
  this->sigma_c         .setDefaultValue(sigma_c);

  this->sigma_c_eff     .reset();
  this->delta_c         .reset();
  this->delta_max       .reset();
  this->insertion_stress.reset();
  this->sigma_c         .reset();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::resizeCohesiveArrays() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::resizeCohesiveArrays();

  this->sigma_c_eff     .resize();
  this->delta_c         .resize();
  this->insertion_stress.resize();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialCohesiveBilinear);


__END_AKANTU__
