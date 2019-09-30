/**
 * @file   material_damage_iterative_viscoelastic.cc
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 *
 * @brief  Specialization of the class material damage to damage only one gauss
 * point at a time and propagate damage in a linear way. Max principal stress
 * criterion is used as a failure criterion. Material is combined with the
 * viscoelastic maxwell material.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_damage_iterative_viscoelastic.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
MaterialDamageIterativeViscoelastic<spatial_dimension>::
    MaterialDamageIterativeViscoelastic(SolidMechanicsModel & model,
                                        const ID & id)
    : parent(model, id), sigma_v_conv("sigma_v_conv", *this),
      epsilon_v_conv("epsilon_v_conv", *this), gradu_last("gradu_last", *this),
      dissipated_energy_damage("dissipated_energy_damage", *this) {

  AKANTU_DEBUG_IN();

  this->dissipated_energy_damage.initialize(1);
  this->dissipated_energy.initializeHistory();
  this->mechanical_work.initializeHistory();
  this->integral.initializeHistory();
  this->damage.initializeHistory();
  this->gradu_last.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  parent::initMaterial();

  UInt stress_size = spatial_dimension * spatial_dimension;
  this->sigma_v.initializeHistory();
  this->epsilon_v.initializeHistory();
  this->sigma_v_conv.initialize(stress_size * this->Ev.size());
  this->epsilon_v_conv.initialize(stress_size * this->Ev.size());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::computeTangentModuli(const ElementType & el_type,
                                             Array<Real> & tangent_matrix,
                                             GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real> & dam = this->damage(el_type);
  auto dam_it = dam.begin();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent, *dam_it);
  ++dam_it;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::afterSolveStep() {

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it =
        this->sigma_v(el_type, _not_ghost)
            .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    auto epsilon_v_it =
        this->epsilon_v(el_type, _not_ghost)
            .begin(spatial_dimension, spatial_dimension, this->Eta.size());
    auto damage_it = this->damage(el_type, _not_ghost).begin();

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    this->updateIntVarOnQuad(grad_u, *previous_gradu_it, *sigma_v_it,
                             *epsilon_v_it, *damage_it);

    ++previous_gradu_it;
    ++sigma_v_it;
    ++epsilon_v_it;
    ++damage_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

    this->updateDissipatedEnergy(el_type);
    this->updateDissipatedEnergyDamage(el_type);
    this->damage.saveCurrentValues();
    this->gradu_last.copy(this->gradu);
  }
}

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::updateIntVariables() {

  this->sigma_v_conv.copy(this->sigma_v);
  this->epsilon_v_conv.copy(this->epsilon_v);
  this->dissipated_energy.saveCurrentValues();
  this->integral.saveCurrentValues();
  this->mechanical_work.saveCurrentValues();
  this->gradu.saveCurrentValues();
  this->stress.saveCurrentValues();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::beforeSolveStep() {

  this->sigma_v.saveCurrentValues();
  this->epsilon_v.saveCurrentValues();
  this->sigma_v.copy(this->sigma_v_conv);
  this->epsilon_v.copy(this->epsilon_v_conv);
  this->integral.copy(this->integral.previous());
  this->dissipated_energy.copy(this->dissipated_energy.previous());
  this->mechanical_work.copy(this->mechanical_work.previous());
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto sigma_v_it =
      this->sigma_v(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto dam_it = this->damage(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeStressOnQuad(grad_u, *previous_gradu_it, sigma, *sigma_v_it,
                      *sigma_th_it, *dam_it);
  ++sigma_th_it;
  ++previous_gradu_it;
  ++sigma_v_it;
  ++dam_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(this->gradu(el_type, ghost_type),
                                          el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->norm_av_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Compute the dissipated energy due to the damage only. For this difference in
 * elastic strain is taken into account. Trapezoidal rule is used. */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::updateDissipatedEnergyDamage(ElementType el_type) {

  auto epsilon_p = this->gradu_last(el_type)
                       .begin(spatial_dimension, spatial_dimension);
  auto sigma_v_it =
      this->sigma_v(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto epsilon_v_it =
      this->epsilon_v(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto sigma_v_pr_it =
      this->sigma_v.previous(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto epsilon_v_pr_it =
      this->epsilon_v.previous(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto dam_it = this->damage(el_type).begin();
  auto dam_pr_it = this->damage.previous(el_type).begin();

  auto epot = this->potential_energy(el_type).begin();
  auto ints = this->int_sigma(el_type).begin();
  auto edd = this->dissipated_energy_damage(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  updateDissipatedEnergyDamageOnQuad(
      grad_u, *epsilon_p, *sigma_v_it, *epsilon_v_it, *sigma_v_pr_it,
      *epsilon_v_pr_it, *dam_it, *dam_pr_it, *epot, *ints, *edd);

  ++sigma_v_it;
  ++epsilon_v_it;
  ++sigma_v_pr_it;
  ++epsilon_v_pr_it;
  ++epsilon_p;
  ++dam_it;
  ++dam_pr_it;
  ++epot;
  ++ints;
  ++edd;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::computePotentialEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computePotentialEnergy(el_type);

  auto epot = this->potential_energy(el_type).begin();
  auto sigma_v_it =
      this->sigma_v(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto epsilon_v_it =
      this->epsilon_v(el_type)
          .begin(spatial_dimension, spatial_dimension, this->Eta.size());
  auto dam_it = this->damage(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  this->computePotentialEnergyOnQuad(grad_u, *epot, *sigma_v_it, *epsilon_v_it,
                                     *dam_it);
  ++epot;
  ++sigma_v_it;
  ++epsilon_v_it;
  ++dam_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialDamageIterativeViscoelastic<spatial_dimension>::getEnergy(
    const std::string & type) {
  if (type == "dissipated_damage")
    return getDissipatedEnergyDamage();
  else
    return MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialDamageIterativeViscoelastic<spatial_dimension>::getEnergy(
    const std::string & energy_id, ElementType type, UInt index) {
  if (energy_id == "dissipated_damage")
    return getDissipatedEnergyDamage(type, index);
  else
    return MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(
        energy_id, type, index);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialDamageIterativeViscoelastic<
    spatial_dimension>::getDissipatedEnergyDamage() const {
  AKANTU_DEBUG_IN();

  Real dde = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    dde += this->fem.integrate(this->dissipated_energy_damage(type, _not_ghost),
                               type, _not_ghost,
                               this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return dde;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialDamageIterativeViscoelastic<
    spatial_dimension>::getDissipatedEnergyDamage(ElementType type,
                                                  UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);
  auto it = this->dissipated_energy_damage(type, _not_ghost)
                .begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[index], type, gindex);
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(damage_iterative_viscoelastic,
                     MaterialDamageIterativeViscoelastic);

} // namespace akantu
