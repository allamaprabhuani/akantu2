/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension, template <UInt> class MatParent>
MaterialVreePeerlings<spatial_dimension, MatParent>::MaterialVreePeerlings(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id), MaterialVreePeerlingsParent(model, id),
      Kapa("Kapa", *this),
      strain_rate_vreepeerlings("strain-rate-vreepeerlings", *this),
      Full_dam_value_strain("fulldam-valstrain", *this),
      Full_dam_value_strain_rate("fulldam-valstrain-rate", *this),
      Number_damage("number-damage", *this), equi_strain("equi-strain", *this),
      equi_strain_rate("equi-strain-rate", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("Kapaoi", Kapaoi, 0.0001, _pat_parsable);
  this->registerParam("Kapac", Kapac, 0.0002, _pat_parsable);
  this->registerParam("Arate", Arate, 0., _pat_parsable);
  this->registerParam("Brate", Brate, 1., _pat_parsable);
  this->registerParam("Crate", Brate, 1., _pat_parsable);
  this->registerParam("Kct", Kct, 1., _pat_parsable);
  this->registerParam("Kapao_randomness", Kapao_randomness, 0., _pat_parsable);

  this->Kapa.initialize(1);
  this->equi_strain.initialize(1);
  this->equi_strain_rate.initialize(1);
  this->Full_dam_value_strain.initialize(1);
  this->Full_dam_value_strain_rate.initialize(1);
  this->Number_damage.initialize(1);
  this->strain_rate_vreepeerlings.initialize(spatial_dimension *
                                             spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlings<spatial_dimension, MatParent>::initMaterial() {
  AKANTU_DEBUG_IN();
  MaterialVreePeerlingsParent::initMaterial();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension, template <UInt> class MatParent>
void MaterialVreePeerlings<spatial_dimension, MatParent>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialVreePeerlingsParent::computeStress(el_type, ghost_type);

  Real * dam = this->damage(el_type, ghost_type).data();
  Real * equi_straint = equi_strain(el_type, ghost_type).data();
  Real * equi_straint_rate = equi_strain_rate(el_type, ghost_type).data();
  Real * Kapaq = Kapa(el_type, ghost_type).data();
  Real * FullDam_Valstrain = Full_dam_value_strain(el_type, ghost_type).data();
  Real * FullDam_Valstrain_rate =
      Full_dam_value_strain_rate(el_type, ghost_type).data();
  Real * Nb_damage = Number_damage(el_type, ghost_type).data();
  Real dt = this->model.getTimeStep();

  Array<Idx> & elem_filter = this->element_filter(el_type, ghost_type);
  Array<Real> & velocity = this->model.getVelocity();
  Array<Real> & strain_rate_vrplgs =
      this->strain_rate_vreepeerlings(el_type, ghost_type);

  this->model.getFEEngine().gradientOnIntegrationPoints(
      velocity, strain_rate_vrplgs, spatial_dimension, el_type, ghost_type,
      elem_filter);

  Array<Real>::matrix_iterator strain_rate_vrplgs_it =
      strain_rate_vrplgs.begin(spatial_dimension, spatial_dimension);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & strain_rate = *strain_rate_vrplgs_it;

  computeStressOnQuad(grad_u, sigma, *dam, *equi_straint, *equi_straint_rate,
                      *Kapaq, dt, strain_rate, *FullDam_Valstrain,
                      *FullDam_Valstrain_rate, *Nb_damage);
  ++dam;
  ++equi_straint;
  ++equi_straint_rate;
  ++Kapaq;
  ++strain_rate_vrplgs_it;
  ++FullDam_Valstrain;
  ++FullDam_Valstrain_rate;
  ++Nb_damage;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
