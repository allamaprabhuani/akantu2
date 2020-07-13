/**
 * @file   material_damage_iterative_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of inline functions of the material damage iterative
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "data_accessor.hh"
#include "material_damage_iterative.hh"
#include "material_elastic_orthotropic_heterogeneous.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_DAMAGE_ITERATIVE_INLINE_IMPL_CC__
#define __AKANTU_MATERIAL_DAMAGE_ITERATIVE_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
MaterialDamageIterative<spatial_dimension, ElasticParent>::
    MaterialDamageIterative(SolidMechanicsModel & model, const ID & id)
    : parent(model, id), Sc("Sc", *this),
      equivalent_stress("equivalent_stress", *this),
      reduction_step("reduction_step", *this),
      crack_normals("crack_normals", *this),
      damage_stored("damage_stored", *this),
      reduction_step_stored("reduction_step_stored", *this),
      extra_volume("extra_volume", *this), max_reductions(0) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam", prescribed_dam, 0.1,
                      _pat_parsable | _pat_modifiable, "prescribed damage");
  this->registerParam(
      "dam_threshold", dam_threshold, 0.8, _pat_parsable | _pat_modifiable,
      "damage threshold at which damage damage will be set to 1");
  this->registerParam(
      "dam_tolerance", dam_tolerance, 0.01, _pat_parsable | _pat_modifiable,
      "damage tolerance to decide if quadrature point will be damaged");
  this->registerParam("max_damage", max_damage, 0.99999,
                      _pat_parsable | _pat_modifiable, "maximum damage value");
  this->registerParam("max_reductions", max_reductions, UInt(10),
                      _pat_parsable | _pat_modifiable, "max reductions");
  this->registerParam("compute_extra_volume", compute_extra_volume, false,
                      _pat_parsmod,
                      "Compute additional volume within cracked elements");

  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->reduction_step.initialize(1);
  this->crack_normals.initialize(spatial_dimension * spatial_dimension);
  this->damage_stored.initialize(1);
  this->reduction_step_stored.initialize(1);
  this->extra_volume.setDefaultValue(0.);
  this->extra_volume.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterative<spatial_dimension, ElasticParent>::
    computeNormalizedEquivalentStress(ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  for (auto && data : zip(make_view(this->stress(el_type, ghost_type),
                                    spatial_dimension, spatial_dimension),
                          make_view(Sc(el_type, ghost_type)),
                          make_view(equivalent_stress(el_type, ghost_type)),
                          make_view(crack_normals(el_type, ghost_type),
                                    spatial_dimension, spatial_dimension))) {

    const auto & sigma = std::get<0>(data);
    const auto & sigma_crit = std::get<1>(data);
    auto & sigma_eq = std::get<2>(data);
    auto & crack_norm = std::get<3>(data);

    /// compute eigenvalues and eigenvectors and sort them
    sigma.eig(eigenvalues, crack_norm, true);

    sigma_eq = eigenvalues[0] / sigma_crit;

    /// normalize each eigenvector
    for (auto && vec : crack_norm) {
      Vector<Real> vector(vec);
      vector.normalize();
    }
    crack_norm = crack_norm.transpose();
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterative<spatial_dimension, ElasticParent>::
    computeAllStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum and average equivalent stresses
  if (ghost_type == _not_ghost) {
    norm_max_equivalent_stress = 0;
  }

  parent::computeAllStresses(ghost_type);

  /// find global Gauss point with highest stress
  // auto rve_model = dynamic_cast<SolidMechanicsModelRVE *>(&this->model);
  // if (rve_model == NULL) {
  /// is no RVE model
  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(norm_max_equivalent_stress, SynchronizerOperation::_max);
  //}
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterative<spatial_dimension, ElasticParent>::
    findMaxNormalizedEquivalentStress(ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (ghost_type == _not_ghost) {

    const Array<Real> & e_stress = equivalent_stress(el_type);
    auto equivalent_stress_it = e_stress.begin();
    auto equivalent_stress_end = e_stress.end();
    Array<Real> & dam = this->damage(el_type);
    auto dam_it = dam.begin();

    for (; equivalent_stress_it != equivalent_stress_end;
         ++equivalent_stress_it, ++dam_it) {
      /// check if max equivalent stress for this element type is greater than
      /// the current norm_max_eq_stress and if the element is not already fully
      /// damaged
      if (*equivalent_stress_it > norm_max_equivalent_stress &&
          *dam_it < max_damage) {
        norm_max_equivalent_stress = *equivalent_stress_it;
      }
    }
  }
  AKANTU_DEBUG_OUT();
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
UInt MaterialDamageIterative<spatial_dimension, ElasticParent>::updateDamage() {
  UInt nb_damaged_elements = 0;
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  if (norm_max_equivalent_stress >= 1.) {

    AKANTU_DEBUG_IN();

    /// update the damage only on non-ghosts elements! Doesn't make sense to
    /// update on ghost.
    auto && element_types =
        this->model.getFEEngine().getMesh().elementTypes(spatial_dimension);
    for (auto el_type : element_types) {
      const auto & equivalent_stresses = equivalent_stress(el_type);
      auto & reduction_steps = this->reduction_step(el_type);
      auto & damages = this->damage(el_type);

      for (auto && data : zip(equivalent_stresses, damages, reduction_steps)) {
        const auto & equivalent_stress = std::get<0>(data);
        auto & dam = std::get<1>(data);
        auto & reduction = std::get<2>(data);

        /// check if damage occurs
        if (equivalent_stress >=
            (1 - dam_tolerance) * norm_max_equivalent_stress) {
          /// check if this element can still be damaged
          if (reduction == this->max_reductions)
            continue;
          reduction += 1;
          if (reduction == this->max_reductions) {
            dam = max_damage;
          } else {
            dam += prescribed_dam;
          }
          nb_damaged_elements += 1;
        }
      }
    }
  }

  // auto * rve_model = dynamic_cast<SolidMechanicsModelRVE
  // *>(&this->model); if (rve_model == NULL) {
  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(nb_damaged_elements, SynchronizerOperation::_sum);
  //}

  AKANTU_DEBUG_OUT();
  return nb_damaged_elements;
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterative<spatial_dimension, ElasticParent>::
    updateEnergiesAfterDamage(ElementType el_type) {
  parent::updateEnergies(el_type);
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
UInt MaterialDamageIterative<spatial_dimension, ElasticParent>::
    updateDamageOnQuad(UInt quad_index, const Real /*eq_stress*/,
                       const ElementType & el_type,
                       const GhostType & ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
                      "Your prescribed damage must be greater than zero");

  Array<Real> & dam = this->damage(el_type, ghost_type);
  Real & dam_on_quad = dam(quad_index);

  /// check if damage occurs
  if (equivalent_stress(el_type, ghost_type)(quad_index) >=
      (1 - dam_tolerance) * norm_max_equivalent_stress) {
    if (dam_on_quad < dam_threshold)
      dam_on_quad += prescribed_dam;
    else
      dam_on_quad = max_damage;
    return 1;
  }

  return 0;
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline UInt
MaterialDamageIterative<spatial_dimension, ElasticParent>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {

  if (tag == SynchronizationTag::_user_2) {
    return sizeof(Real) * this->getModel().getNbIntegrationPoints(elements);
  }

  return parent::getNbData(elements, tag);
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline void MaterialDamageIterative<spatial_dimension, ElasticParent>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_user_2) {
    DataAccessor<Element>::packElementalDataHelper(
        this->damage, buffer, elements, true, this->damage.getFEEngine());
  }

  return parent::packData(buffer, elements, tag);
}

/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline void
MaterialDamageIterative<spatial_dimension, ElasticParent>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_user_2) {
    DataAccessor<Element>::unpackElementalDataHelper(
        this->damage, buffer, elements, true, this->damage.getFEEngine());
  }
  return parent::unpackData(buffer, elements, tag);
}
} // namespace akantu

/* --------------------------------------------------------------------------
 */

#endif /* __AKANTU_MATERIAL_DAMAGE_ITERATIVE_INLINE_IMPL_CC__ */
