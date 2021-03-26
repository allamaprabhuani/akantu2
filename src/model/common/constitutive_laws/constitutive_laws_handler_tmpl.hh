/**
 * @file   constitutive_laws_handler_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation  ven mar 26 2021
 *
 * @brief Implementation of the consistutive laws handler.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "constitutive_law_non_local_interface.hh"
#include "constitutive_laws_handler.hh"
#include "non_local_manager.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH
#define AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline decltype(auto)
ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaws() {
  return make_transform_adaptor(
      constitutive_laws, [](auto && value) -> decltype(auto) {
        return aka::as_type<ConstitutiveLawType>(*value);
      });
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline decltype(auto)
ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaws() const {
  return make_transform_adaptor(
      constitutive_laws, [](auto && value) -> decltype(auto) {
        return aka::as_type<ConstitutiveLawType>(*value);
      });
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline UInt
ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLawIndex(
    const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << parent_id << " has no constitutive_law named " << name);
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline auto & ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaw(
    UInt cl_index) {
  AKANTU_DEBUG_ASSERT(cl_index < constitutive_laws.size(),
                      "The model " << parent_id
                                   << " has no constitutive law number "
                                   << cl_index);
  return aka::as_type<ConstitutiveLawType>(*constitutive_laws.at(cl_index));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline const auto &
ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaw(
    UInt cl_index) const {
  AKANTU_DEBUG_ASSERT(cl_index < constitutive_laws.size(),
                      "The model " << parent_id
                                   << " has no constitutive law number "
                                   << cl_index);
  return aka::as_type<ConstitutiveLawType>(*constitutive_laws.at(cl_index));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline auto & ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaw(
    const std::string & name) {
  std::map<std::string, UInt>::const_iterator it =
      constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << parent_id << " has no constitutive_law named " << name);
  }

  return aka::as_type<ConstitutiveLawType>(*constitutive_laws[it->second]);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
inline const auto &
ConstitutiveLawsHandler<ConstitutiveLawType>::getConstitutiveLaw(
    const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << parent_id << " has no constitutive_law named " << name);
  }
  return aka::as_type<ConstitutiveLawType>(*constitutive_laws[it->second]);
}

/* -------------------------------------------------------------------------- */
#define AKA_FWD_CL(...) ::std::forward<decltype(__VA_ARGS__)>(__VA_ARGS__)

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
template <typename Operation>
void ConstitutiveLawsHandler<ConstitutiveLawType>::splitByConstitutiveLaw(
    const Array<Element> & elements, Operation && op) const {
  std::vector<Array<Element>> elements_per_cl(constitutive_laws.size());
  this->splitElementByConstitutiveLaw(elements, elements_per_cl);

  for (auto && cl : zip(constitutive_laws, elements_per_cl)) {
    AKA_FWD_CL(op)(AKA_FWD_CL(*std::get<0>(cl)), AKA_FWD_CL(std::get<1>(cl)));
  }
}

#undef AKA_FWD_CL

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
auto & ConstitutiveLawsHandler<ConstitutiveLawType>::registerNewConstitutiveLaw(
    const ParserSection & section) {
  std::string cl_name;
  std::string cl_type = section.getName();
  std::string opt_param = section.getOption();

  try {
    std::string tmp = section.getParameter("name");
    cl_name = tmp; /** this can seam weird, but there is an ambiguous operator
                    * overload that i couldn't solve. @todo remove the
                    * weirdness of this code
                    */
  } catch (debug::Exception &) {
    AKANTU_ERROR("A constitutive_law of type \'"
                 << cl_type
                 << "\' in the input file has been defined without a name!");
  }

  auto && cl = this->registerNewConstitutiveLaw(cl_name, cl_type, opt_param);

  cl.parseSection(section);

  return cl;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
auto & ConstitutiveLawsHandler<ConstitutiveLawType>::registerNewConstitutiveLaw(
    const ID & cl_name, const ID & cl_type, const ID & opt_param) {
  AKANTU_DEBUG_ASSERT(constitutive_laws_names_to_id.find(cl_name) ==
                          constitutive_laws_names_to_id.end(),
                      "A constitutive_law with this name '"
                          << cl_name << "' has already been registered. "
                          << "Please use unique names for constitutive_laws");

  UInt cl_count = constitutive_laws.size();
  constitutive_laws_names_to_id[cl_name] = cl_count;

  std::stringstream sstr_cl;
  sstr_cl << this->id << ":" << cl_count << ":" << cl_type;
  ID cl_id = sstr_cl.str();

  std::unique_ptr<ConstitutiveLaw> constitutive_law =
      ConstitutiveLawType::getFactory().allocate(cl_type, _spatial_dimension,
                                                 opt_param, *this, cl_id);

  constitutive_laws.push_back(std::move(constitutive_law));

  return aka::as_type<ConstitutiveLawType>(*(constitutive_laws.back()));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<
    ConstitutiveLawType>::instantiateConstitutiveLaws() {
  ParserSection model_section;
  bool is_empty;
  std::tie(model_section, is_empty) = this->getParserSection();

  if (not is_empty) {
    auto model_constitutive_laws =
        model_section.getSubSections(ParserType::_constitutive_law);
    for (const auto & section : model_constitutive_laws) {
      this->registerNewConstitutiveLaw(section);
    }
  }

  auto sub_sections =
      this->parser.getSubSections(ParserType::_constitutive_law);
  for (const auto & section : sub_sections) {
    this->registerNewConstitutiveLaw(section);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  for (auto & constitutive_law : constitutive_laws) {
    if (dynamic_cast<ConstitutiveLawNonLocalInterface *>(
            constitutive_law.get()) == nullptr) {
      continue;
    }

    this->non_local_manager = std::make_unique<NonLocalManager>(
        *this, *this, parent_id + ":non_local_manager");
    break;
  }
#endif

  if (constitutive_laws.empty()) {
    AKANTU_EXCEPTION("No constitutive_laws where instantiated for the model"
                     << parent_id);
  }
  are_constitutive_laws_instantiated = true;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::
    assignConstitutiveLawToElements(const ElementTypeMapArray<UInt> * filter) {

  for_each_element(
      _mesh,
      [&](auto && element) {
        UInt cl_index = (*constitutive_law_selector)(element);
        AKANTU_DEBUG_ASSERT(cl_index < constitutive_laws.size(),
                            "The constitutive_law selector returned an index "
                            "that does not exists");
        constitutive_law_index(element) = cl_index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  if (non_local_manager) {
    non_local_manager->synchronize(*this,
                                   SynchronizationTag::_constitutive_law_id);
  }

  for_each_element(
      _mesh,
      [&](auto && element) {
        auto cl_index = constitutive_law_index(element);
        auto index = constitutive_laws[cl_index]->addElement(element);
        constitutive_law_local_numbering(element) = index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  // synchronize the element constitutive_law arrays
  this->synchronize(SynchronizationTag::_constitutive_law_id);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::initConstitutiveLaws() {
  AKANTU_DEBUG_ASSERT(not constitutive_laws.empty(),
                      "No constitutive_law to initialize !");

  this->assignConstitutiveLawToElements();

  for (auto & constitutive_law : constitutive_laws) {
    /// init internals properties
    constitutive_law->initConstitutiveLaw();
  }

  this->synchronize(SynchronizationTag::_smm_init_cl);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
Int ConstitutiveLawsHandler<ConstitutiveLawType>::getInternalIndexFromID(
    const ID & id) const {
  AKANTU_DEBUG_IN();

  auto it = constitutive_laws.begin();
  auto end = constitutive_laws.end();

  for (; it != end; ++it) {
    if ((*it)->getID() == id) {
      AKANTU_DEBUG_OUT();
      return (it - constitutive_laws.begin());
    }
  }

  AKANTU_DEBUG_OUT();
  return -1;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::reassignConstitutiveLaw() {
  AKANTU_DEBUG_IN();

  std::vector<Array<Element>> element_to_add(constitutive_laws.size());
  std::vector<Array<Element>> element_to_remove(constitutive_laws.size());

  for_each_element(mesh, [&](auto && element) {
    auto old_constitutive_law = constitutive_law_index(element);
    auto new_constitutive_law = (*constitutive_law_selector)(element);
    if (old_constitutive_law != new_constitutive_law) {
      element_to_add[new_constitutive_law].push_back(element);
      element_to_remove[old_constitutive_law].push_back(element);
    }
  });

  for (auto && data : enumerate(constitutive_laws)) {
    auto cl_index = std::get<0>(data);
    auto & cl = *std::get<1>(data);

    cl.removeElements(element_to_remove[cl_index]);
    cl.addElements(element_to_add[cl_index]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::
    updateDataForNonLocalCriterion(ElementTypeMapReal & criterion) {
  const ID field_name = criterion.getName();
  for (auto & constitutive_law : constitutive_laws) {
    if (!constitutive_law->isInternal<Real>(field_name, _ek_regular)) {
      continue;
    }

    for (auto ghost_type : ghost_types) {
      constitutive_law->flattenInternal(field_name, criterion, ghost_type,
                                        _ek_regular);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::initializeNonLocal() {
  this->non_local_manager->synchronize(
      *this, SynchronizationTag::_constitutive_law_id);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::
    insertIntegrationPointsInNeighborhoods(GhostType ghost_type) {
  for (auto & cl : constitutive_laws) {
    ConstitutiveLawNonLocalInterface * cl_non_local;
    if ((cl_non_local = dynamic_cast<ConstitutiveLawNonLocalInterface *>(
             cl.get())) == nullptr) {
      continue;
    }

    ElementTypeMapArray<Real> quadrature_points_coordinates(
        "quadrature_points_coordinates_tmp_nl", this->id);
    quadrature_points_coordinates.initialize(this->getFEEngine(),
                                             _nb_component = spatial_dimension,
                                             _ghost_type = ghost_type);

    for (const auto & type : quadrature_points_coordinates.elementTypes(
             Model::spatial_dimension, ghost_type)) {
      this->getFEEngine().computeIntegrationPointsCoordinates(
          quadrature_points_coordinates(type, ghost_type), type, ghost_type);
    }

    cl_non_local->initConstitutiveLawNonLocal();

    cl_non_local->insertIntegrationPointsInNeighborhoods(
        ghost_type, quadrature_points_coordinates);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::updateLocalInternal(
    ElementTypeMapReal & internal_flat, GhostType ghost_type,
    ElementKind kind) {
  const ID field_name = internal_flat.getName();
  for (auto & constitutive_law : constitutive_laws) {
    if (constitutive_law->isInternal<Real>(field_name, kind)) {
      constitutive_law->flattenInternal(field_name, internal_flat, ghost_type,
                                        kind);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::updateNonLocalInternal(
    ElementTypeMapReal & internal_flat, GhostType ghost_type,
    ElementKind kind) {

  const ID field_name = internal_flat.getName();

  for (auto & mat : constitutive_laws) {
    if (not aka::is_of_type<ConstitutiveLawNonLocalInterface>(*mat)) {
      continue;
    }

    auto & mat_non_local =
        dynamic_cast<ConstitutiveLawNonLocalInterface &>(*mat);
    mat_non_local.updateNonLocalInternals(internal_flat, field_name, ghost_type,
                                          kind);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::
    splitElementByConstitutiveLaw(
        const Array<Element> & elements,
        std::vector<Array<Element>> & elements_per_cl) const {
  for (const auto & el : elements) {
    Element cl_el = el;
    cl_el.element = this->constitutive_law_local_numbering(el);
    elements_per_cl[this->constitutive_law_index(el)].push_back(cl_el);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
UInt ConstitutiveLawsHandler<ConstitutiveLawType>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
    UInt size{0};

    if(tag == SynchronizationTag::_constitutive_law_id) {
        size += elements.size() * sizeof(UInt);
        break;
    }


  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      size += cl.getNbData(elements, tag);
    });
  }
  return size;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_constitutive_law_id) {
    DataAccessor<Element>::packElementalDataHelper(
        constitutive_law_index, buffer, elements, false, getFEEngine());
  }

  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      cl.packData(buffer, elements, tag);
    });
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_constitutive_law_id) {
    for (auto && element : elements) {
      UInt recv_cl_index;
      buffer >> recv_cl_index;
      UInt & cl_index = constitutive_law_index(element);
      if (cl_index != UInt(-1)) {
        continue;
      }

      // add ghosts element to the correct constitutive_law
      cl_index = recv_cl_index;
      UInt index = constitutive_laws[cl_index]->addElement(element);
      constitutive_law_local_numbering(element) = index;
    }
  }

  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      cl.unpackData(buffer, elements, tag);
    });
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  this->constitutive_law_index.initialize(mesh, _element_kind = _ek_not_defined,
                                          _with_nb_element = true,
                                          _default_value = UInt(-1));
  this->constitutive_law_local_numbering.initialize(
      mesh, _element_kind = _ek_not_defined, _with_nb_element = true,
      _default_value = UInt(-1));

  ElementTypeMapArray<UInt> filter("new_element_filter", this->getID());

  for (const auto & elem : element_list) {
    if (mesh.getSpatialDimension(elem.type) != spatial_dimension) {
      continue;
    }

    if (!filter.exists(elem.type, elem.ghost_type)) {
      filter.alloc(0, 1, elem.type, elem.ghost_type);
    }
    filter(elem.type, elem.ghost_type).push_back(elem.element);
  }

  // this fails in parallel if the event is sent on facet between constructor
  // and initFull \todo: to debug...
  this->assignConstitutiveLawToElements(&filter);

  for (auto & constitutive_law : constitutive_laws) {
    constitutive_law->onElementsAdded(element_list, event);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType>
void ConstitutiveLawsHandler<ConstitutiveLawType>::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<UInt> & new_numbering,
    const RemovedElementsEvent & event) {
  for (auto & constitutive_law : constitutive_laws) {
    constitutive_law->onElementsRemoved(element_list, new_numbering, event);
  }
}

} // namespace akantu

#endif /* AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH */
