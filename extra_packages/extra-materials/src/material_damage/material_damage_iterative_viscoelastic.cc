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
    : iterative_parent(model, id), sigma_v_conv("sigma_v_conv", *this),
      epsilon_v_conv("epsilon_v_conv", *this), gradu_last("gradu_last", *this),
      dissipated_energy_damage("dissipated_energy_damage", *this) {

  AKANTU_DEBUG_IN();

  // this->dissipated_energy_damage.initialize(1);
  // this->dissipated_energy.initializeHistory();
  // this->mechanical_work.initializeHistory();
  // this->integral.initializeHistory();
  // this->damage.initializeHistory();
  // this->gradu_last.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  iterative_parent::initMaterial();

  UInt stress_size = spatial_dimension * spatial_dimension;
  this->sigma_v.initializeHistory();
  // this->epsilon_v.initializeHistory();
  this->sigma_v_conv.initialize(stress_size * this->Ev.size());
  // this->epsilon_v_conv.initialize(stress_size * this->Ev.size());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::computeTangentModuli(ElementType el_type,
                                             Array<Real> & tangent_matrix,
                                             GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real> & dam = this->damage(el_type);
  auto dam_it = dam.begin();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  viscous_grandparent::computeTangentModuliOnQuad(tangent, *dam_it);
  ++dam_it;
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<
    spatial_dimension>::updateIntVariables() {

  this->sigma_v_conv.copy(this->sigma_v);
  // this->epsilon_v_conv.copy(this->epsilon_v);
  // this->dissipated_energy.saveCurrentValues();
  // this->integral.saveCurrentValues();
  // this->mechanical_work.saveCurrentValues();
  this->gradu.saveCurrentValues();
  this->stress.saveCurrentValues();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::beforeSolveStep() {

  this->sigma_v.saveCurrentValues();
  // this->epsilon_v.saveCurrentValues();
  this->sigma_v.copy(this->sigma_v_conv);
  // this->epsilon_v.copy(this->epsilon_v_conv);
  // this->integral.copy(this->integral.previous());
  // this->dissipated_energy.copy(this->dissipated_energy.previous());
  // this->mechanical_work.copy(this->mechanical_work.previous());
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

  viscous_grandparent::computeStressOnQuad(grad_u, *previous_gradu_it, sigma,
                                          *sigma_v_it, *sigma_th_it, *dam_it);
  ++sigma_th_it;
  ++previous_gradu_it;
  ++sigma_v_it;
  ++dam_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  this->computeNormalizedEquivalentStress(el_type, ghost_type);
  this->norm_max_equivalent_stress = 0;
  this->findMaxNormalizedEquivalentStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeViscoelastic<spatial_dimension>::afterSolveStep(bool converged) {

  if (not converged) {
    return;
  }

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it =
        this->sigma_v(el_type, _not_ghost)
            .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    // auto epsilon_v_it =
    //     this->epsilon_v(el_type, _not_ghost)
    //         .begin(spatial_dimension, spatial_dimension, this->Eta.size());
    auto damage_it = this->damage(el_type, _not_ghost).begin();

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    this->updateSigmaViscOnQuad(grad_u, *previous_gradu_it, *sigma_v_it,
                                /**epsilon_v_it,*/ *damage_it);

    ++previous_gradu_it;
    ++sigma_v_it;
    // ++epsilon_v_it;
    ++damage_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

    // this->updateDissipatedEnergy(el_type);
    // this->updateDissipatedEnergyDamage(el_type);
    // this->damage.saveCurrentValues();
    // this->gradu_last.copy(this->gradu);
  }
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(damage_iterative_viscoelastic,
                     MaterialDamageIterativeViscoelastic);

} // namespace akantu
