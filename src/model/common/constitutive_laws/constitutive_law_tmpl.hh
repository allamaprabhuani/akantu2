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
#include "fe_engine.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAW_TMPL_HH
#define AKANTU_CONSTITUTIVE_LAW_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T = Real,
          template <typename T> InternalFieldType = InternalField>
inline std::shared_ptr<InternalFieldType<T>>
ConstitutiveLawInternalHandler::registerInternal<T>(const ID & id,
                                                    UInt nb_component) {
  auto && internal = std::make_shared<InternalFieldType<T>>(
      id, *this, dim, this->fe_engine_id, this->element_filter);
  internal->initialize(nb_component);
  internal_vectors[internal->getRegisterID()] = internal;
  return internal;
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::unregisterInternal(const ID & id) {
  internal_vectors.erase(id);
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::savePreviousState() {
  for (auto pair : internal_vectors) {
    if (pair.second->hasHistory()) {
      pair.second->saveCurrentValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::restorePreviousState() {
  for (auto pair : internal_vectors) {
    if (pair.second->hasHistory()) {
      pair.second->restorePreviousValues();
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::resizeInternals() {
  for (auto && pair : internal_vectors) {
    pair.second->resize();
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
const InternalField<T> &
ConstitutiveLawInternalHandler::getInternal(const ID & id) const {
  auto it = internal_vectors.find(this->id + ":" + id);
  if (it != internal_vectors.end() and
      aka::is_of_type<InternalField<T>>(*it->second)) {
    return aka::as_type<InternalField<T>>(*it->second);
  }

  AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                          << ") does not contain an internal "
                                          << id << " (" << (getID() + ":" + id)
                                          << ")");
}

/* -------------------------------------------------------------------------- */
template <typename T>
InternalField<T> & ConstitutiveLawInternalHandler::getInternal(const ID & id) {
  auto it = internal_vectors.find(getID() + ":" + id);
  if (it != internal_vectors.end() and
      aka::is_of_type<InternalField<T>>(*it->second)) {
    return aka::as_type<InternalField<T>>(*it->second);
  }

  AKANTU_SILENT_EXCEPTION("The material " << name << "(" << getID()
                                          << ") does not contain an internal "
                                          << id << " (" << (getID() + ":" + id)
                                          << ")");
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool ConstitutiveLawInternalHandler::isInternal(
    const ID & id, const ElementKind & element_kind) const {
  auto it = internal_vectors.find(this->getID() + ":" + id);

  return (it != internal_vectors.end() and
          aka::is_of_type<InternalField<T>>(*it->second) and
          aka::as_type<InternalField<T>>(*it->second).getElementKind() !=
              element_kind);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename T>
const Array<T> &
ConstitutiveLawInternalHandler::getArray(const ID & vect_id, ElementType type,
                                         GhostType ghost_type) const {
  if (isInternal<T>(vect_id, Mesh::getKind(type))) {
    auto && internal = this->template getInternal<T>(vect_id);
    if (internal.exists(type, ghost_type)) {
      return internal(type, ghost_type);
    }

    AKANTU_SILENT_EXCEPTION(
        "The internal " << vect_id << " in the constitutive law " << name
                        << " (" << getID() << ") does not contain the type ["
                        << type << ":" << ghost_type << "]");
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << " (" << getID()
                          << ") does not contain an internal field "
                          << vect_id);
}

/* -------------------------------------------------------------------------- */
template <typename T>
Array<T> & ConstitutiveLawInternalHandler::getArray(const ID & vect_id,
                                                    ElementType type,
                                                    GhostType ghost_type) {
  if (isInternal<T>(vect_id, Mesh::getKind(type))) {
    auto && internal = this->template getInternal<T>(vect_id);
    if (internal.exists(type, ghost_type)) {
      return internal(type, ghost_type);
    }
    AKANTU_SILENT_EXCEPTION(
        "The internal " << vect_id << " in the constitutive law " << name
                        << " (" << getID() << ") does not contain the type ["
                        << type << ":" << ghost_type << "]");
  }

  AKANTU_SILENT_EXCEPTION("The constitutive law "
                          << name << " (" << getID()
                          << ") does not contain an internal field "
                          << vect_id);
}

/* -------------------------------------------------------------------------- */
inline void ConstitutiveLawInternalHandler::removeIntegrationPoints(
    ElementTypeMapArray<UInt> & new_numbering) {
  for (auto pair : internal_vectors) {
    pair.second->removeIntegrationPoints(new_numbering);
  }
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
ConstitutiveLaw<ConstitutiveLawsHandler_>::ConstitutiveLaw(
    ConstitutiveLawsHandler_ & handler, const ID & id, UInt spatial_dimension,
    ElementKind element_kind, const ID & fe_engine_id)
    : ConstitutiveLawInternalHandler(id, spatial_dimension),
      Parsable(handler.getConstitutiveLawParserType(), id), handler(handler),
      element_filter("element_filter", id), default_fe_engine_id(fe_engine_id) {

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

  this->updateInternalParameters();

  is_init = true;
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

  this->removeIntegrationPoints(constitutive_law_local_new_numbering);

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

  this->removeIntegrationPoints(constitutive_law_local_new_numbering);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void ConstitutiveLaw<ConstitutiveLawsHandler_>::packInternalFieldHelper(
    const InternalField<T> & data_to_pack, CommunicationBuffer & buffer,
    const Array<Element> & elements) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements,
                                           data_to_pack.getFEEngine());
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
inline void
ConstitutiveLaw<ConstitutiveLawsHandler_>::unpackInternalFieldHelper(
    InternalField<T> & data_to_unpack, CommunicationBuffer & buffer,
    const Array<Element> & elements) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements,
                                             data_to_unpack.getFEEngine());
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
inline Element ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToLocalElement(
    const Element & global_element) const {
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
ConstitutiveLaw<ConstitutiveLawsHandler_>::convertToGlobalElement(
    const Element & local_element) const {
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
template <typename T>
inline ElementTypeMap<UInt>
ConstitutiveLaw<ConstitutiveLawsHandler_>::getInternalDataPerElem(
    const ID & field_id, ElementKind element_kind) const {

  if (not this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field "
                     << id << " in the constitutive law " << this->name);
  }

  const auto & internal_field = this->template getInternal<T>(field_id);
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
template <typename T>
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

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawsHandler_>
template <typename T>
void ConstitutiveLaw<ConstitutiveLawsHandler_>::inflateInternal(
    const std::string & field_id, const ElementTypeMapArray<T> & field,
    GhostType ghost_type, ElementKind element_kind) {
  if (!this->template isInternal<T>(field_id, element_kind)) {
    AKANTU_EXCEPTION("Cannot find internal field " << id << " in material "
                                                   << this->name);
  }

  InternalField<T> & internal_field = this->template getInternal<T>(field_id);
  const FEEngine & fe_engine = internal_field.getFEEngine();

  for (auto && type : field.elementTypes(ghost_type)) {
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

} // namespace akantu

#endif // AKANTU_CONSTITUTIVE_LAW_TMPL_HH
