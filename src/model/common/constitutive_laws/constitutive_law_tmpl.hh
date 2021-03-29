/**
 * @file   constitutive_law_tmpl.hh *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Implementation of the templated part of the constitutive law class
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
#include "constitutive_law.hh"
#include "constitutive_laws_handler.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAW_TMPL_HH
#define AKANTU_CONSTITUTIVE_LAW_TMPL_HH

namespace akantu {
template <class ConstitutiveLawsHandler_>
ConstitutiveLaw<ConstitutiveLawsHandler_>::ConstitutiveLaw(
    ConstitutiveLawsHandler_ & handler, const ID & id, ElementKind element_kind)
    : Parsable(ParserType::_constitutive_law, id), handler(handler), id(id),
      element_filter("element_filter", id) {
  /// for each connectivity types allocate the element filer array of the
  /// constitutive law
  element_filter.initialize(handler.getMesh(),
                            _spatial_dimension = handler.getSpatialDimension(),
                            _element_kind = element_kind);

  this->initialize();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::initialize() {
  registerParam("name", name, std::string(), _pat_parsable | _pat_readable);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::initConstitutiveLaw() {
  this->resizeInternals();

  this->updateInternalParmaters();

  is_init = true;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::savePreviousState() {
  for (auto pair : internal_vectors_real) {
    if (pair.second->hasHistory()) {
      pair.second->saveCurrentValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::restorePreviousState() {
  for (auto pair : internal_vectors_real) {
    if (pair.second->hasHistory()) {
      pair.second->restorePreviousValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::resizeInternals() {
  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it) {
    it->second->resize();
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->resize();
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::addElements(
    const Array<Element> & elements_to_add) {
  AKANTU_DEBUG_IN();

  UInt law_id = handler.getConstitutiveLawIndex(name);
  for (const auto & element : elements_to_add) {
    auto index = this->addElement(element);
    handler.constitutive_law_index(element) = law_id;
    handler.constitutive_law_local_numbering(element) = index;
  }

  this->resizeInternals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::removeElements(
    const Array<Element> & elements_to_remove) {
  AKANTU_DEBUG_IN();

  auto el_begin = elements_to_remove.begin();
  auto el_end = elements_to_remove.end();

  if (elements_to_remove.empty()) {
    return;
  }

  auto & mesh = handler.getMesh();

  ElementTypeMapArray<UInt> constitutive_law_local_new_numbering(
      "remove constitutive law filter elem", id);

  constitutive_law_local_new_numbering.initialize(
      mesh, _element_filter = &element_filter, _element_kind = _ek_not_defined,
      _with_nb_element = true);

  ElementTypeMapArray<UInt> element_filter_tmp("element_filter_tmp", id);

  element_filter_tmp.initialize(mesh, _element_filter = &element_filter,
                                _element_kind = _ek_not_defined);

  ElementTypeMap<UInt> new_ids, element_ids;

  for_each_element(
      mesh,
      [&](auto && el) {
        if (not new_ids(el.type, el.ghost_type)) {
          element_ids(el.type, el.ghost_type) = 0;
        }

        auto & element_id = element_ids(el.type, el.ghost_type);
        auto l_el = Element{el.type, element_id, el.ghost_type};
        if (std::find(el_begin, el_end, el) != el_end) {
          constitutive_law_local_new_numbering(l_el) = UInt(-1);
          return;
        }

        element_filter_tmp(el.type, el.ghost_type).push_back(el.element);
        if (not new_ids(el.type, el.ghost_type)) {
          new_ids(el.type, el.ghost_type) = 0;
        }

        auto & new_id = new_ids(el.type, el.ghost_type);

        constitutive_law_local_new_numbering(l_el) = new_id;
        handler.constitutive_law_local_numbering(el) = new_id;

        ++new_id;
        ++element_id;
      },
      _element_filter = &element_filter, _element_kind = _ek_not_defined);

  for (auto ghost_type : ghost_types) {
    for (const auto & type : element_filter.elementTypes(
             _ghost_type = ghost_type, _element_kind = _ek_not_defined)) {
      element_filter(type, ghost_type)
          .copy(element_filter_tmp(type, ghost_type));
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::onElementsAdded(
    const Array<Element> & /*unused*/, const NewElementsEvent & /*unused*/) {
  this->resizeInternals();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & /*event*/) {

  UInt my_num = handler.getInternalIndexFromID(getID());

  ElementTypeMapArray<UInt> constitutive_law_local_new_numbering(
      "remove constitutive law filter elem", getID());

  auto el_begin = element_list.begin();
  auto el_end = element_list.end();

  for (auto && gt : ghost_types) {
    for (auto && type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {

      if (not element_filter.exists(type, gt) ||
          element_filter(type, gt).empty()) {
        continue;
      }

      auto & elem_filter = element_filter(type, gt);
      auto & law_indexes = this->handler.constitutive_law_index(type, gt);
      auto & law_loc_num =
          this->handler.constitutive_law_local_numbering(type, gt);
      auto nb_element = this->handler.getMesh().getNbElement(type, gt);

      // all constitutive laws will resized to the same size...
      law_indexes.resize(nb_element);
      law_loc_num.resize(nb_element);

      if (!constitutive_law_local_new_numbering.exists(type, gt)) {
        constitutive_law_local_new_numbering.alloc(elem_filter.size(), 1, type,
                                                   gt);
      }

      auto & law_renumbering = constitutive_law_local_new_numbering(type, gt);
      const auto & renumbering = new_numbering(type, gt);
      Array<UInt> elem_filter_tmp;
      UInt ni = 0;
      Element el{type, 0, gt};

      for (UInt i = 0; i < elem_filter.size(); ++i) {
        el.element = elem_filter(i);
        if (std::find(el_begin, el_end, el) == el_end) {
          UInt new_el = renumbering(el.element);
          AKANTU_DEBUG_ASSERT(new_el != UInt(-1),
                              "A not removed element as been badly renumbered");
          elem_filter_tmp.push_back(new_el);
          law_renumbering(i) = ni;

          law_indexes(new_el) = my_num;
          law_loc_num(new_el) = ni;
          ++ni;
        } else {
          law_renumbering(i) = UInt(-1);
        }
      }

      elem_filter.resize(elem_filter_tmp.size());
      elem_filter.copy(elem_filter_tmp);
    }
  }

  for (auto it = internal_vectors_real.begin();
       it != internal_vectors_real.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }

  for (auto it = internal_vectors_uint.begin();
       it != internal_vectors_uint.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }

  for (auto it = internal_vectors_bool.begin();
       it != internal_vectors_bool.end(); ++it) {
    it->second->removeIntegrationPoints(constitutive_law_local_new_numbering);
  }
}


/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void ConstitutiveLaw<ConstitutiveLawsHandler_>::packElementDataHelper(
    const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements, true,
                                           handler.getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void ConstitutiveLaw<ConstitutiveLawsHandler_>::unpackElementDataHelper(
    ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements, const ID & fem_id) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             true, handler.getFEEngine(fem_id));
}  


/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Element
ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToLocalElement(const Element & global_element) const {
  UInt ge = global_element.element;
#ifndef AKANTU_NDEBUG
  UInt model_law_index = handler.getConstitutiveLawByElement(
      global_element.type, global_element.ghost_type)(ge);

  UInt law_index = handler.getConstitutiveLawIndex(this->name);
  AKANTU_DEBUG_ASSERT(model_law_index == law_index,
                      "Conversion of a global  element in a local element for "
                      "the wrong constitutive law "
                          << this->name << std::endl);
#endif
  UInt le = handler.getConstitutiveLawLocalNumbering(
      global_element.type, global_element.ghost_type)(ge);

  Element tmp_quad{global_element.type, le, global_element.ghost_type};
  return tmp_quad;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Element
ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToGlobalElement(const Element & local_element) const {
  UInt le = local_element.element;
  UInt ge =
      this->element_filter(local_element.type, local_element.ghost_type)(le);

  Element tmp_quad{local_element.type, ge, local_element.ghost_type};
  return tmp_quad;
}
  
  

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline UInt ConstitutiveLaw<ConstitutiveLawsHandler_>::addElement(
    ElementType type, UInt element, GhostType ghost_type) {
  Array<UInt> & el_filter = this->element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline UInt
ConstitutiveLaw<ConstitutiveLawsHandler_>::addElement(const Element & element) {
  return this->addElement(element.type, element.element, element.ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline const Parameter &
ConstitutiveLaw<ConstitutiveLawsHandler_>::getParam(const ID & param) const {
  try {
    return get(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the constitutive law"
                                     << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::setParam(const ID & param, T value) {
  try {
    set<T>(param, value);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the constitutive law "
                                     << getID());
  }
  updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::registerInternal<Real, fps>(
    InternalField<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::registerInternal<UInt, fps>(
    InternalField<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}

template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::registerInternal<bool, fps>(
    InternalField<bool> & vect) {
  internal_vectors_bool[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::unregisterInternal<Real, fps>(
    InternalField<Real> & vect) {
  internal_vectors_real.erase(vect.getID());
}

template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::unregisterInternal<UInt, fps>(
    InternalField<UInt> & vect) {
  internal_vectors_uint.erase(vect.getID());
}

template <class ConstitutiveLawsHandler_>
template <int fps>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::unregisterInternal<bool, fps>(
    InternalField<bool> & vect) {
  internal_vectors_bool.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
const InternalField<Real> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<Real, fps>(
    const ID & int_id) const {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <>
InternalField<Real, fps> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<Real, fps>(
    const ID & int_id) {
  auto it = internal_vectors_real.find(getID() + ":" + int_id);
  if (it == internal_vectors_real.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
const InternalField<UInt> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<UInt, fps>(
    const ID & int_id) const {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
InternalField<UInt> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<UInt, fps>(
    const ID & int_id) {
  auto it = internal_vectors_uint.find(getID() + ":" + int_id);
  if (it == internal_vectors_uint.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
const InternalField<bool> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<bool, fps>(
    const ID & int_id) const {
  auto it = internal_vectors_bool.find(getID() + ":" + int_id);
  if (it == internal_vectors_bool.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
InternalField<bool> &
ConstituitveLaw<ConstitutiveLawsHandler_>::getInternal<bool, fps>(
    const ID & int_id) {
  auto it = internal_vectors_bool.find(getID() + ":" + int_id);
  if (it == internal_vectors_bool.end()) {
    AKANTU_SILENT_EXCEPTION("The constitutive law " << name << "(" << getID()
                                            << ") does not contain an internal "
                                            << int_id << " ("
                                            << (getID() + ":" + int_id) << ")");
  }
  return *it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
inline bool ConstitutiveLaw<ConstitutiveLawsHandler_>::isInternal<Real, fps>(
    const ID & id, ElementKind element_kind) const {
  auto internal_array = internal_vectors_real.find(this->getID() + ":" + id);

  return not(internal_array == internal_vectors_real.end() ||
             internal_array->second->getElementKind() != element_kind);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
inline bool ConstitutiveLaw<ConstitutiveLawsHandler_>::isInternal<UInt, fps>(
    const ID & id, ElementKind element_kind) const {
  auto internal_array = internal_vectors_uint.find(this->getID() + ":" + id);

  return not(internal_array == internal_vectors_uint.end() ||
             internal_array->second->getElementKind() != element_kind);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <int fps>
inline bool ConstitutiveLaw<ConstitutiveLawsHandler_>::isInternal<bool, fps>(
    const ID & id, ElementKind element_kind) const {
  auto internal_array = internal_vectors_bool.find(this->getID() + ":" + id);

  return not(internal_array == internal_vectors_bool.end() ||
             internal_array->second->getElementKind() != element_kind);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline ElementTypeMap<UInt>
ConstitutiveLaw<ConstitutiveLawsHandler_>::getInternalDataPerElem(
    const ID & field_id, ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field "
                     << id << " in the constitutive law " << this->name);
  }

  const InternalField<T> & internal_field =
      this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();
  UInt nb_data_per_quad = internal_field.getNbComponent();

  ElementTypeMap<UInt> res;
  for (auto ghost_type : ghost_types) {
    for (auto & type : internal_field.elementTypes(ghost_type)) {
      UInt nb_quadrature_points =
          fe_engine.getNbIntegrationPoints(type, ghost_type);
      res(type, ghost_type) = nb_data_per_quad * nb_quadrature_points;
    }
  }

  return res;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::flattenInternal(
    const std::string & field_id, ElementTypeMapArray<T> & internal_flat,
    const GhostType ghost_type, ElementKind element_kind) const {

  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field "
                     << id << " in the constitutive law " << this->name);
  }

  const InternalField<T> & internal_field =
      this->template getInternal<T>(field_id);

  const FEEngine & fe_engine = internal_field.getFEEngine();
  const Mesh & mesh = fe_engine.getMesh();

  for (auto && type : internal_field.filterTypes(ghost_type)) {
    const Array<Real> & src_vect = internal_field(type, ghost_type);
    const Array<UInt> & filter = internal_field.getFilter(type, ghost_type);

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

    if (nb_element_src == 0) {
      continue;
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


} // namespace akantu

#endif // AKANTU_CONSTITUTIVE_LAW_TMPL_HH
