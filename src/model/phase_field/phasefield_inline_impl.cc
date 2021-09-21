/**
 * @file   phasefield_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 19 2020
 * @date last modification: Fri Apr 02 2021
 *
 * @brief  Phase field implementation of inline functions
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_INLINE_IMPL_CC__
#define __AKANTU_PHASEFIELD_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::addElement(const ElementType & type, UInt element,
                                   const GhostType & ghost_type) {
  Array<UInt> & el_filter = this->element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::addElement(const Element & element) {
  return this->addElement(element.type, element.element, element.ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
PhaseField::registerInternal<Real>(InternalPhaseField<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template <>
inline void
PhaseField::registerInternal<UInt>(InternalPhaseField<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}

template <>
inline void
PhaseField::registerInternal<bool>(InternalPhaseField<bool> & vect) {
  internal_vectors_bool[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template <>
inline void
PhaseField::unregisterInternal<Real>(InternalPhaseField<Real> & vect) {
  internal_vectors_real.erase(vect.getID());
}

template <>
inline void
PhaseField::unregisterInternal<UInt>(InternalPhaseField<UInt> & vect) {
  internal_vectors_uint.erase(vect.getID());
}

template <>
inline void
PhaseField::unregisterInternal<bool>(InternalPhaseField<bool> & vect) {
  internal_vectors_bool.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool PhaseField::isInternal(__attribute__((unused)) const ID & id,
                                   __attribute__((unused))
                                   const ElementKind & element_kind) const {
  AKANTU_TO_IMPLEMENT();
}

template <>
inline bool
PhaseField::isInternal<Real>(const ID & id,
                             const ElementKind & element_kind) const {
  auto internal_array = internal_vectors_real.find(this->getID() + ":" + id);

  return !(internal_array == internal_vectors_real.end() ||
           internal_array->second->getElementKind() != element_kind);
}

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::getNbData(__attribute__((unused))
                                  const Array<Element> & elements,
                                  __attribute__((unused))
                                  const SynchronizationTag & tag) const {
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void PhaseField::packData(__attribute__((unused))
                                 CommunicationBuffer & buffer,
                                 __attribute__((unused))
                                 const Array<Element> & elements,
                                 __attribute__((unused))
                                 const SynchronizationTag & tag) const {}

/* -------------------------------------------------------------------------- */
inline void
PhaseField::unpackData(__attribute__((unused)) CommunicationBuffer & buffer,
                       __attribute__((unused)) const Array<Element> & elements,
                       __attribute__((unused)) const SynchronizationTag & tag) {
}

/* -------------------------------------------------------------------------- */
inline const Parameter & PhaseField::getParam(const ID & param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material "
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void PhaseField::packElementDataHelper(
    const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements, true,
                                           model.getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void PhaseField::unpackElementDataHelper(
    ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             true, model.getFEEngine(fem_id));
}

} // namespace akantu

#endif
