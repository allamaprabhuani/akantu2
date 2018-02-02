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
template<UInt spatial_dimension>
inline void
MaterialDamageIterative<spatial_dimension>::computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam) {
  sigma *= 1-dam;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
UInt MaterialDamageIterative<spatial_dimension>::updateDamage(UInt quad_index, const Real eq_stress, const ElementType & el_type, const GhostType & ghost_type) {
  AKANTU_DEBUG_ASSERT(prescribed_dam > 0.,
		      "Your prescribed damage must be greater than zero");


  Array<Real> & dam = this->damage(el_type, ghost_type);
  Real & dam_on_quad   = dam(quad_index);

  /// check if damage occurs
  if (equivalent_stress(el_type, ghost_type)(quad_index) >= (1-dam_tolerance) * norm_max_equivalent_stress) {
    if (dam_on_quad < dam_threshold)
      dam_on_quad +=prescribed_dam;
    else dam_on_quad = max_damage;
    return 1;
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt MaterialDamageIterative<spatial_dimension>::getNbDataForElements(const Array<Element> & elements,
									     SynchronizationTag tag) const {

  if (tag == _gst_user_2) {
    return sizeof(Real) * this->getModel().getNbIntegrationPoints(elements);
  }

  return MaterialDamage<spatial_dimension>::getNbDataForElements(elements, tag);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialDamageIterative<spatial_dimension>::packElementData(CommunicationBuffer & buffer,
									const Array<Element> & elements,
									SynchronizationTag tag) const {
  if (tag == _gst_user_2) {
    this->packElementDataHelper(this->damage, buffer, elements);
  }

  return MaterialDamage<spatial_dimension>::packElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialDamageIterative<spatial_dimension>::unpackElementData(CommunicationBuffer & buffer,
									  const Array<Element> & elements,
									  SynchronizationTag tag) {
  if (tag == _gst_user_2) {
    this->unpackElementDataHelper(this->damage, buffer, elements);
  }

  return MaterialDamage<spatial_dimension>::unpackElementData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
