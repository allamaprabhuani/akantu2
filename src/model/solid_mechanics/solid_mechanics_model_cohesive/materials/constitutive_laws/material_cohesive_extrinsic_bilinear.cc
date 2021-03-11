/**
 * @file   material_cohesive_linear.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Shenghan Zhang <shenghan.zhang@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Dec 1 2016
 * @date last modification: Fri Mar 5 2021
 *
 * @brief  Bilinear Cohesive  law for extrinsic elements
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_cohesive_extrinsic_bilinear.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveExtrinsicBilinear<spatial_dimension>::MaterialCohesiveExtrinsicBilinear(
    SolidMechanicsModel & model, const ID & id)
  : MaterialCohesiveLinear<spatial_dimension>(model, id),
  delta_h_eff("delta_h_eff", *this) {

  AKANTU_DEBUG_IN();
  
  this->registerParam("h", h, Real(0.),
		      _pat_parsmod, "Slope of the inflection point");
  this->registerParam("delta_h", delta_h, Real(0.),
                      _pat_parsable | _pat_readable, "Inflection displacement");

  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExtrinsicBilinear<spatial_dimension>::initMaterial() {
  MaterialCohesiveLinear<spatial_dimension>::initMaterial();
  this->delta_h_eff.initialize(1);
  this->delta_h_eff.setDefaultValue(delta_h);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveExtrinsicBilinear<spatial_dimension>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {

  /// define iterators
  auto traction_it = this->tractions(el_type, ghost_type).begin(spatial_dimension);
  auto traction_end = this->tractions(el_type, ghost_type).end(spatial_dimension);
  auto opening_it = this->opening(el_type, ghost_type).begin(spatial_dimension);
  auto contact_traction_it =
    this->contact_tractions(el_type, ghost_type).begin(spatial_dimension);
  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(spatial_dimension);

  auto normal_it = normal.begin(spatial_dimension);
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto delta_max_prev_it =
    this->delta_max.previous(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto delta_h_it = this->delta_h_eff(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto insertion_stress_it =
    this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it,
         ++sigma_c_it, ++delta_max_it, ++delta_c_it, ++damage_it,
         ++insertion_stress_it, ++delta_max_prev_it, ++delta_h_it,
	 ++contact_opening_it, ++contact_traction_it) {

    Real normal_opening_norm{0};
    Real tangential_opening_norm{0};
    bool penetration{false};
    
    this->computeTractionOnQuad(
        *traction_it, *opening_it, *normal_it, *delta_max_it, *delta_c_it,
	*delta_h_it, *insertion_stress_it, *sigma_c_it,
	normal_opening, tangential_opening, normal_opening_norm,
	tangential_opening_norm, *damage_it, penetration,
        *contact_traction_it, *contact_opening_it);
    
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(cohesive_extrinsic_bilinear, MaterialCohesiveExtrinsicBilinear);

} // akantu
