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

#ifndef __AKANTU_PHASEFIELD_INLINE_IMPL_HH__
#define __AKANTU_PHASEFIELD_INLINE_IMPL_HH__

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
template <typename T>
void PhaseField::flattenInternal(const std::string & field_id,
                                 ElementTypeMapArray<T> & internal_flat,
                                 const GhostType ghost_type,
                                 ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in phasefield "
                                                   << this->name);
  }

  const InternalPhaseField<T> & internal_field =
      this->template getInternal<T>(field_id);

  const FEEngine & fe_engine = internal_field.getFEEngine();
  const Mesh & mesh = fe_engine.getMesh();

  for (auto && type : internal_field.filterTypes(ghost_type)) {
    const auto & src_vect = internal_field(type, ghost_type);
    const auto & filter = internal_field.getFilter(type, ghost_type);

    // total number of elements in the corresponding mesh
    UInt nb_element_dst = mesh.getNbElement(type, ghost_type);
    // number of element in the internal field
    UInt nb_element_src = filter.size();
    // number of quadrature points per elem
    UInt nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    // number of data per quadrature point
    UInt nb_data_per_quad = internal_field.getNbComponent();

    if (!internal_flat.exists(type, ghost_type)) {
      internal_flat.alloc(nb_element_dst * nb_quad_per_elem, nb_data_per_quad,
                          type, ghost_type);
    }

    // number of data per element
    UInt nb_data = nb_quad_per_elem * nb_data_per_quad;

    Array<Real> & dst_vect = internal_flat(type, ghost_type);
    dst_vect.resize(nb_element_dst * nb_quad_per_elem);

    auto it_dst = make_view(dst_vect, nb_data).begin();

    for (auto && data : zip(filter, make_view(src_vect, nb_data))) {
      it_dst[std::get<0>(data)] = std::get<1>(data);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
void PhaseField::inflateInternal(const std::string & field_id,
                                 const ElementTypeMapArray<T> & field,
                                 GhostType ghost_type,
                                 ElementKind element_kind) {
  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in phasefield "
                                                   << this->name);
  }

  InternalPhaseField<T> & internal_field =
      this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();

  for (auto && type : field.elementTypes(spatial_dimension, ghost_type)) {
    if (not internal_field.exists(type, ghost_type)) {
      continue;
    }
    const auto & filter = internal_field.getFilter(type, ghost_type);

    const auto & src_array = field(type, ghost_type);
    auto & dest_array = internal_field(type, ghost_type);

    auto nb_quad_per_elem = fe_engine.getNbIntegrationPoints(type);
    auto nb_component = src_array.getNbComponent();

    AKANTU_DEBUG_ASSERT(
        field.size() == fe_engine.getMesh().getNbElement(type, ghost_type) *
                            nb_quad_per_elem,
        "The ElementTypeMapArray to inflate is not of the proper size");
    AKANTU_DEBUG_ASSERT(
        dest_array.getNbComponent() == nb_component,
        "The ElementTypeMapArray has not the proper number of components");

    auto src =
        make_view(field(type, ghost_type), nb_component, nb_quad_per_elem)
            .begin();
    for (auto && data :
         zip(filter, make_view(dest_array, nb_component, nb_quad_per_elem))) {
      std::get<1>(data) = src[std::get<0>(data)];
    }
  }
}

/* -------------------------------------------------------------------------- */
inline UInt PhaseField::getNbData(const Array<Element> & elements,
                                  const SynchronizationTag & tag) const {

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void PhaseField::packData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) const {}

/* -------------------------------------------------------------------------- */
inline void PhaseField::unpackData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag & tag) {}

/* -------------------------------------------------------------------------- */
inline const Parameter & PhaseField::getParam(const ID & param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the phasefield "
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

/* -------------------------------------------------------------------------- */
template <typename T>
inline const InternalPhaseField<T> &
PhaseField::getInternal([[gnu::unused]] const ID & int_id) const {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline InternalPhaseField<T> &
PhaseField::getInternal([[gnu::unused]] const ID & int_id) {
  AKANTU_TO_IMPLEMENT();
  return NULL;
}

/* -------------------------------------------------------------------------- */
template <>
inline const InternalPhaseField<Real> &
PhaseField::getInternal(const ID & int_id) const {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The phasefield "
                            << name << "(" << getID()
                            << ") does not contain an internal " << int_id
                            << " (" << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline InternalPhaseField<Real> & PhaseField::getInternal(const ID & int_id) {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The phasefield "
                            << name << "(" << getID()
                            << ") does not contain an internal " << int_id
                            << " (" << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline const InternalPhaseField<UInt> &
PhaseField::getInternal(const ID & int_id) const {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The phasefield "
                            << name << "(" << getID()
                            << ") does not contain an internal " << int_id
                            << " (" << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <>
inline InternalPhaseField<UInt> & PhaseField::getInternal(const ID & int_id) {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The phasefield "
                            << name << "(" << getID()
                            << ") does not contain an internal " << int_id
                            << " (" << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const Array<T> & PhaseField::getArray(const ID & vect_id,
                                             ElementType type,
                                             GhostType ghost_type) const {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The phasefield " << name << "(" << getID()
                                              << ") does not contain a vector "
                                              << vect_id << " [" << e << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> & PhaseField::getArray(const ID & vect_id, ElementType type,
                                       GhostType ghost_type) {
  try {
    return this->template getInternal<T>(vect_id)(type, ghost_type);
  } catch (debug::Exception & e) {
    AKANTU_SILENT_EXCEPTION("The phasefield " << name << "(" << getID()
                                              << ") does not contain a vector "
                                              << vect_id << " [" << e << "]");
  }
}

} // namespace akantu

#endif
