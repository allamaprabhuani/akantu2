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
      reduction_step("reduction_step", *this),
      equivalent_stress("equivalent_stress", *this), max_reductions(0),
      min_equivalent_stress("min_equivalent_stress", *this),
      crack_normals("normal_to_crack", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("Sc", Sc, _pat_parsable, "critical stress threshold");
  this->registerParam("prescribed_dam", prescribed_dam, 0.1,
                      _pat_parsable | _pat_modifiable, "prescribed damage");
  this->registerParam(
      "dam_threshold", dam_threshold, 0.8, _pat_parsable | _pat_modifiable,
      "damage threshold at which damage damage will be set to 1");
  this->registerParam(
      "dam_tolerance", dam_tolerance, 0.01, _pat_parsable | _pat_modifiable,
      "damage tolerance to decide if quadrature point will be damageed");
  this->registerParam("max_damage", max_damage, 0.99999,
                      _pat_parsable | _pat_modifiable, "maximum damage value");
  this->registerParam("max_reductions", max_reductions, UInt(10),
                      _pat_parsable | _pat_modifiable, "max reductions");
  this->registerParam("contact", contact, false,
                      _pat_parsable | _pat_modifiable,
                      "parameter responsible for the stiffness recovery when "
                      "in compression");
  this->registerParam("smoothen_stiffness_change", smoothen_stiffness_change,
                      false, _pat_parsable | _pat_modifiable,
                      "smoothening stiffness change from damaged to "
                      "non-damaged when in compression");
  this->registerParam("K", K, 5000., _pat_parsable | _pat_modifiable,
                      "parameter of smoothening law");

  this->use_previous_stress = true;
  this->use_previous_gradu = true;
  this->Sc.initialize(1);
  this->equivalent_stress.initialize(1);
  this->reduction_step.initialize(1);
  this->min_equivalent_stress.initialize(1);
  this->crack_normals.initialize(spatial_dimension * spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
void MaterialDamageIterative<spatial_dimension, ElasticParent>::
    computeNormalizedEquivalentStress(const Array<Real> & /*grad_us*/,
                                      ElementType el_type,
                                      GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// Vector to store eigenvalues of current stress tensor
  Vector<Real> eigenvalues(spatial_dimension);

  for (auto && data : zip(make_view(this->stress(el_type, ghost_type),
                                    spatial_dimension, spatial_dimension),
                          make_view(Sc(el_type, ghost_type)),
                          make_view(equivalent_stress(el_type, ghost_type)),
                          make_view(min_equivalent_stress(el_type, ghost_type)),
                          make_view(crack_normals(el_type, ghost_type),
                                    spatial_dimension, spatial_dimension))) {

    const auto & sigma = std::get<0>(data);
    const auto & sigma_crit = std::get<1>(data);
    auto & sigma_eq = std::get<2>(data);
    auto & min_sigma_eq = std::get<3>(data);
    auto & crack_norm = std::get<4>(data);

    /// compute eigenvalues and eigenvectors and sort them
    sigma.eig(eigenvalues, crack_norm, true);

    sigma_eq = eigenvalues[0] / sigma_crit;
    min_sigma_eq = eigenvalues[spatial_dimension - 1] / sigma_crit;

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
/* -----------------------------------------------------------------------*/
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline auto MaterialDamageIterative<spatial_dimension, ElasticParent>::
    computePrincStrainAndRotMatrix(const Matrix<Real> & sigma,
                                   const Matrix<Real> & grad_u,
                                   bool max_strain) {

  Vector<Real> eigEps(spatial_dimension);
  Matrix<Real> rotation_matrix(spatial_dimension, spatial_dimension);

  // small strain tensor
  Matrix<Real> strain(grad_u);
  strain += grad_u.transpose();
  strain *= 0.5;

  // compute eigenvalues
  strain.eig(eigEps, rotation_matrix, false);

  /// normalize each column of the rotation matrix by the length of
  /// corresponding eigen vector
  for (auto && c : arange(rotation_matrix.cols())) {
    Vector<Real> vect(rotation_matrix(c));
    vect /= vect.norm();
  }

  Matrix<Real> RtSigma(spatial_dimension, spatial_dimension);
  Matrix<Real> SigmaPrime(spatial_dimension, spatial_dimension);
  Vector<Real> eigSigma(spatial_dimension);
  RtSigma.mul<true, false>(rotation_matrix, sigma);
  SigmaPrime.mul<false, false>(RtSigma, rotation_matrix);
  for (auto i : arange(spatial_dimension))
    eigSigma[i] = SigmaPrime(i, i);

  // position of the biggest or smallest stress value
  Real pos;
  if (max_strain)
    pos =
        std::distance(eigSigma.storage(),
                      std::max_element(eigSigma.storage(),
                                       eigSigma.storage() + spatial_dimension));
  else
    pos =
        std::distance(eigSigma.storage(),
                      std::min_element(eigSigma.storage(),
                                       eigSigma.storage() + spatial_dimension));

  return std::make_tuple(eigEps, eigSigma, rotation_matrix, pos);
}
/* --------------------------------------------------------------------------
 */
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline Real MaterialDamageIterative<spatial_dimension, ElasticParent>::
    computeSmoothingFactor(const Real & eps, const Real & sigma_prime,
                           const Real & dam, const Real & delta0) {
  // smoothening between two slopes is done by tanh function
  // as the center of smoothening minus delta0 is taken
  Real smooth_coef = 1. - dam;
  if (sigma_prime < 0)
    smooth_coef = 1. - dam / 2 * (1 + tanh(this->K * (eps + delta0)));
  return smooth_coef;
}
/* ----------------------------------------------------------------------*/
template <UInt spatial_dimension, template <UInt> class ElasticParent>
inline void
MaterialDamageIterative<spatial_dimension, ElasticParent>::rotateTensor(
    Matrix<Real> & T, const Matrix<Real> & rotation_matrix) {
  AKANTU_DEBUG_ASSERT(T.rows() == rotation_matrix.rows() and
                          T.cols() == rotation_matrix.cols(),
                      "Dimensions of tensors do not match");

  Matrix<Real> temp(T);
  temp.mul<true, false>(rotation_matrix, T);
  T.mul<false, false>(temp, rotation_matrix);
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
