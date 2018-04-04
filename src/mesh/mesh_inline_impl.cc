/**
 * @file   mesh_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Implementation of the inline functions of the mesh class
 *
 * @section LICENSE
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
#include "aka_iterators.hh"
#include "element_class.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_INLINE_IMPL_CC__
#define __AKANTU_MESH_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
inline ElementKind Element::kind() const { return Mesh::getKind(type); }

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <typename... pack>
Mesh::ElementTypesIteratorHelper Mesh::elementTypes(pack &&... _pack) const {
  return connectivities.elementTypes(_pack...);
}

/* -------------------------------------------------------------------------- */
inline RemovedNodesEvent::RemovedNodesEvent(const Mesh & mesh)
    : new_numbering(mesh.getNbNodes(), 1, "new_numbering") {}

/* -------------------------------------------------------------------------- */
inline RemovedElementsEvent::RemovedElementsEvent(const Mesh & mesh,
                                                  const ID & new_numbering_id)
    : new_numbering(new_numbering_id, mesh.getID(), mesh.getMemoryID()) {}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<NewElementsEvent>(NewElementsEvent & event) {
  this->nodes_to_elements.resize(nodes->size());
  for (const auto & elem : event.getList()) {
    const Array<UInt> & conn = connectivities(elem.type, elem.ghost_type);

    UInt nb_nodes_per_elem = this->getNbNodesPerElement(elem.type);

    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      UInt node = conn(elem.element, n);
      if (not nodes_to_elements[node])
        nodes_to_elements[node] = std::make_unique<std::set<Element>>();
      nodes_to_elements[node]->insert(elem);
    }
  }

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <> inline void Mesh::sendEvent<NewNodesEvent>(NewNodesEvent & event) {
  this->computeBoundingBox();

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <>
inline void
Mesh::sendEvent<RemovedElementsEvent>(RemovedElementsEvent & event) {
  this->connectivities.onElementsRemoved(event.getNewNumbering());
  this->fillNodesToElements();
  this->computeBoundingBox();

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<RemovedNodesEvent>(RemovedNodesEvent & event) {
  const auto & new_numbering = event.getNewNumbering();
  this->removeNodesFromArray(*nodes, new_numbering);
  if (nodes_global_ids and not mesh_parent)
    this->removeNodesFromArray(*nodes_global_ids, new_numbering);
  if (nodes_type and not mesh_parent)
    this->removeNodesFromArray(*nodes_type, new_numbering);

  if (not nodes_to_elements.empty()) {
    std::vector<std::unique_ptr<std::set<Element>>> tmp(
        nodes_to_elements.size());
    auto it = nodes_to_elements.begin();

    UInt new_nb_nodes = 0;
    for (auto new_i : new_numbering) {
      if (new_i != UInt(-1)) {
        tmp[new_i] = std::move(*it);
        ++new_nb_nodes;
      }
      ++it;
    }

    tmp.resize(new_nb_nodes);
    std::move(tmp.begin(), tmp.end(), nodes_to_elements.begin());
  }

  computeBoundingBox();

  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Mesh::removeNodesFromArray(Array<T> & vect,
                                       const Array<UInt> & new_numbering) {
  Array<T> tmp(vect.size(), vect.getNbComponent());
  UInt nb_component = vect.getNbComponent();
  UInt new_nb_nodes = 0;
  for (UInt i = 0; i < new_numbering.size(); ++i) {
    UInt new_i = new_numbering(i);
    if (new_i != UInt(-1)) {
      T * to_copy = vect.storage() + i * nb_component;
      std::uninitialized_copy(to_copy, to_copy + nb_component,
                              tmp.storage() + new_i * nb_component);
      ++new_nb_nodes;
    }
  }

  tmp.resize(new_nb_nodes);
  vect.copy(tmp);
}

/* -------------------------------------------------------------------------- */
inline Array<UInt> & Mesh::getNodesGlobalIdsPointer() {
  AKANTU_DEBUG_IN();
  if (not nodes_global_ids) {
    nodes_global_ids = std::make_shared<Array<UInt>>(
        nodes->size(), 1, getID() + ":nodes_global_ids");

    for (auto && global_ids : enumerate(*nodes_global_ids)) {
      std::get<1>(global_ids) = std::get<0>(global_ids);
    }
  }

  AKANTU_DEBUG_OUT();
  return *nodes_global_ids;
}

/* -------------------------------------------------------------------------- */
inline Array<NodeType> & Mesh::getNodesTypePointer() {
  AKANTU_DEBUG_IN();
  if (not nodes_type) {
    nodes_type =
        std::make_shared<Array<NodeType>>(nodes->size(), 1, _nt_normal);
  }

  AKANTU_DEBUG_OUT();
  return *nodes_type;
}

/* -------------------------------------------------------------------------- */
inline Array<UInt> &
Mesh::getConnectivityPointer(const ElementType & type,
                             const GhostType & ghost_type) {
  if (connectivities.exists(type, ghost_type))
    return connectivities(type, ghost_type);

  if (ghost_type != _not_ghost) {
    ghosts_counters.alloc(0, 1, type, ghost_type, 1);
  }

  AKANTU_DEBUG_INFO("The connectivity vector for the type " << type
                                                            << " created");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  return connectivities.alloc(0, nb_nodes_per_element, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline Array<std::vector<Element>> &
Mesh::getElementToSubelementPointer(const ElementType & type,
                                    const GhostType & ghost_type) {
  return getDataPointer<std::vector<Element>>("element_to_subelement", type,
                                              ghost_type, 1, true);
}

/* -------------------------------------------------------------------------- */
inline Array<Element> &
Mesh::getSubelementToElementPointer(const ElementType & type,
                                    const GhostType & ghost_type) {
  auto & array = getDataPointer<Element>(
      "subelement_to_element", type, ghost_type, getNbFacetsPerElement(type),
      true, is_mesh_facets, ElementNull);
  return array;
}

/* -------------------------------------------------------------------------- */
inline const auto & Mesh::getElementToSubelement() const {
  return getData<std::vector<Element>>("element_to_subelement");
}

/* -------------------------------------------------------------------------- */
inline const auto &
Mesh::getElementToSubelement(const ElementType & type,
                             const GhostType & ghost_type) const {
  return getData<std::vector<Element>>("element_to_subelement", type,
                                       ghost_type);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getElementToSubelement(const ElementType & type,
                                           const GhostType & ghost_type) {
  return getData<std::vector<Element>>("element_to_subelement", type,
                                       ghost_type);
}

/* -------------------------------------------------------------------------- */
inline const auto &
Mesh::getElementToSubelement(const Element & element) const {
  return getData<std::vector<Element>>("element_to_subelement")(element);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getElementToSubelement(const Element & element) {
  return getData<std::vector<Element>>("element_to_subelement")(element);
}

/* -------------------------------------------------------------------------- */
inline const auto & Mesh::getSubelementToElement() const {
  return getData<Element>("subelement_to_element");
}

/* -------------------------------------------------------------------------- */
inline const auto &
Mesh::getSubelementToElement(const ElementType & type,
                             const GhostType & ghost_type) const {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getSubelementToElement(const ElementType & type,
                                           const GhostType & ghost_type) {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline VectorProxy<Element>
Mesh::getSubelementToElement(const Element & element) const {
  const auto & sub_to_element =
      this->getSubelementToElement(element.type, element.ghost_type);
  auto it = sub_to_element.begin(sub_to_element.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
inline VectorProxy<Element>
Mesh::getSubelementToElement(const Element & element) {
  auto & sub_to_element =
      this->getSubelementToElement(element.type, element.ghost_type);
  auto it = sub_to_element.begin(sub_to_element.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> &
Mesh::getDataPointer(const ID & data_name, const ElementType & el_type,
                     const GhostType & ghost_type, UInt nb_component,
                     bool size_to_nb_element, bool resize_with_parent) {
  Array<T> & tmp = mesh_data.getElementalDataArrayAlloc<T>(
      data_name, el_type, ghost_type, nb_component);

  if (size_to_nb_element) {
    if (resize_with_parent)
      tmp.resize(mesh_parent->getNbElement(el_type, ghost_type));
    else
      tmp.resize(this->getNbElement(el_type, ghost_type));
  } else {
    tmp.resize(0);
  }

  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> &
Mesh::getDataPointer(const ID & data_name, const ElementType & el_type,
                     const GhostType & ghost_type, UInt nb_component,
                     bool size_to_nb_element, bool resize_with_parent,
                     const T & defaul_) {
  Array<T> & tmp = mesh_data.getElementalDataArrayAlloc<T>(
      data_name, el_type, ghost_type, nb_component);

  if (size_to_nb_element) {
    if (resize_with_parent)
      tmp.resize(mesh_parent->getNbElement(el_type, ghost_type), defaul_);
    else
      tmp.resize(this->getNbElement(el_type, ghost_type), defaul_);
  } else {
    tmp.resize(0);
  }

  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline bool Mesh::hasData(const ID & data_name, const ElementType & el_type,
                          const GhostType & ghost_type) const {
  return mesh_data.hasDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const Array<T> & Mesh::getData(const ID & data_name,
                                      const ElementType & el_type,
                                      const GhostType & ghost_type) const {
  return mesh_data.getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline Array<T> & Mesh::getData(const ID & data_name,
                                const ElementType & el_type,
                                const GhostType & ghost_type) {
  return mesh_data.getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::hasData(const ID & data_name) const {
  return mesh_data.hasData(data_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const ElementTypeMapArray<T> &
Mesh::getData(const ID & data_name) const {
  return mesh_data.getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline ElementTypeMapArray<T> & Mesh::getData(const ID & data_name) {
  return mesh_data.getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline ElementTypeMapArray<T> & Mesh::registerData(const ID & data_name) {
  this->mesh_data.registerElementalData<T>(data_name);
  return this->getData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const ElementType & type,
                               const GhostType & ghost_type) const {
  try {

    const Array<UInt> & conn = connectivities(type, ghost_type);
    return conn.size();
  } catch (...) {
    return 0;
  }
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbElement(const UInt spatial_dimension,
                               const GhostType & ghost_type,
                               const ElementKind & kind) const {
  AKANTU_DEBUG_ASSERT(spatial_dimension <= 3 || spatial_dimension == UInt(-1),
                      "spatial_dimension is " << spatial_dimension
                                              << " and is greater than 3 !");
  UInt nb_element = 0;

  for (auto type : elementTypes(spatial_dimension, ghost_type, kind))
    nb_element += getNbElement(type, ghost_type);

  return nb_element;
}

/* -------------------------------------------------------------------------- */
inline void Mesh::getBarycenter(const Element & element,
                                Vector<Real> & barycenter) const {
  Vector<UInt> conn = getConnectivity(element);
  Matrix<Real> local_coord(spatial_dimension, conn.size());
  auto node_begin = make_view(*nodes, spatial_dimension).begin();

  for(auto && node : enumerate(conn)) {
    local_coord(std::get<0>(node)) = Vector<Real>(node_begin[std::get<1>(node)]);
  }

  Math::barycenter(local_coord.storage(), conn.size(), spatial_dimension,
                   barycenter.storage());
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElement(const ElementType & type) {
  UInt nb_nodes_per_element = 0;
#define GET_NB_NODES_PER_ELEMENT(type)                                         \
  nb_nodes_per_element = ElementClass<type>::getNbNodesPerElement()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_NODES_PER_ELEMENT);
#undef GET_NB_NODES_PER_ELEMENT
  return nb_nodes_per_element;
}

/* -------------------------------------------------------------------------- */
inline ElementType Mesh::getP1ElementType(const ElementType & type) {
  ElementType p1_type = _not_defined;
#define GET_P1_TYPE(type) p1_type = ElementClass<type>::getP1ElementType()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_P1_TYPE);
#undef GET_P1_TYPE
  return p1_type;
}

/* -------------------------------------------------------------------------- */
inline ElementKind Mesh::getKind(const ElementType & type) {
  ElementKind kind = _ek_not_defined;
#define GET_KIND(type) kind = ElementClass<type>::getKind()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_KIND);
#undef GET_KIND
  return kind;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getSpatialDimension(const ElementType & type) {
  UInt spatial_dimension = 0;
#define GET_SPATIAL_DIMENSION(type)                                            \
  spatial_dimension = ElementClass<type>::getSpatialDimension()
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SPATIAL_DIMENSION);
#undef GET_SPATIAL_DIMENSION

  return spatial_dimension;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetTypes(const ElementType & type,
                                  __attribute__((unused)) UInt t) {
  UInt nb = 0;
#define GET_NB_FACET_TYPE(type) nb = ElementClass<type>::getNbFacetTypes()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_FACET_TYPE);
#undef GET_NB_FACET_TYPE
  return nb;
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getFacetType(const ElementType & type, UInt t) {
#define GET_FACET_TYPE(type) return ElementClass<type>::getFacetType(t);

  AKANTU_BOOST_ALL_ELEMENT_SWITCH_NO_DEFAULT(GET_FACET_TYPE);

#undef GET_FACET_TYPE

  return _not_defined;
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getAllFacetTypes(const ElementType & type) {
#define GET_FACET_TYPE(type) return ElementClass<type>::getFacetTypes();

  AKANTU_BOOST_ALL_ELEMENT_SWITCH_NO_DEFAULT(GET_FACET_TYPE);
#undef GET_FACET_TYPE

  return ElementClass<_not_defined>::getFacetTypes();
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElement(const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt n_facet = 0;
#define GET_NB_FACET(type) n_facet = ElementClass<type>::getNbFacetsPerElement()

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_FACET);
#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbFacetsPerElement(const ElementType & type, UInt t) {
  AKANTU_DEBUG_IN();

  UInt n_facet = 0;
#define GET_NB_FACET(type)                                                     \
  n_facet = ElementClass<type>::getNbFacetsPerElement(t)

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NB_FACET);
#undef GET_NB_FACET

  AKANTU_DEBUG_OUT();
  return n_facet;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getFacetLocalConnectivity(const ElementType & type, UInt t) {
  AKANTU_DEBUG_IN();

#define GET_FACET_CON(type)                                                    \
  AKANTU_DEBUG_OUT();                                                          \
  return ElementClass<type>::getFacetLocalConnectivityPerElement(t)

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_FACET_CON);
#undef GET_FACET_CON

  AKANTU_DEBUG_OUT();
  return ElementClass<_not_defined>::getFacetLocalConnectivityPerElement(0);
  // This avoid a compilation warning but will certainly
  // also cause a segfault if reached
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getFacetConnectivity(const Element & element, UInt t) const {
  AKANTU_DEBUG_IN();

  Matrix<const UInt> local_facets(getFacetLocalConnectivity(element.type, t));
  Matrix<UInt> facets(local_facets.rows(), local_facets.cols());

  const Array<UInt> & conn = connectivities(element.type, element.ghost_type);

  for (UInt f = 0; f < facets.rows(); ++f) {
    for (UInt n = 0; n < facets.cols(); ++n) {
      facets(f, n) = conn(element.element, local_facets(f, n));
    }
  }

  AKANTU_DEBUG_OUT();
  return facets;
}

/* -------------------------------------------------------------------------- */
inline VectorProxy<UInt> Mesh::getConnectivity(const Element & element) const {
  const auto & conn = connectivities(element.type, element.ghost_type);
  auto it = conn.begin(conn.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
inline VectorProxy<UInt> Mesh::getConnectivity(const Element & element) {
  auto & conn = connectivities(element.type, element.ghost_type);
  auto it = conn.begin(conn.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Mesh::extractNodalValuesFromElement(
    const Array<T> & nodal_values, T * local_coord, UInt * connectivity,
    UInt n_nodes, UInt nb_degree_of_freedom) const {
  for (UInt n = 0; n < n_nodes; ++n) {
    memcpy(local_coord + n * nb_degree_of_freedom,
           nodal_values.storage() + connectivity[n] * nb_degree_of_freedom,
           nb_degree_of_freedom * sizeof(T));
  }
}

/* -------------------------------------------------------------------------- */
inline void Mesh::addConnectivityType(const ElementType & type,
                                      const GhostType & ghost_type) {
  getConnectivityPointer(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isPureGhostNode(UInt n) const {
  return nodes_type ? ((*nodes_type)(n) == _nt_pure_ghost) : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalOrMasterNode(UInt n) const {
  return nodes_type
             ? ((*nodes_type)(n) == _nt_master) ||
                   ((*nodes_type)(n) == _nt_normal)
             : true;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isLocalNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == _nt_normal : true;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isMasterNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) == _nt_master : false;
}

/* -------------------------------------------------------------------------- */
inline bool Mesh::isSlaveNode(UInt n) const {
  return nodes_type ? (*nodes_type)(n) >= 0 : false;
}

/* -------------------------------------------------------------------------- */
inline NodeType Mesh::getNodeType(UInt local_id) const {
  return nodes_type ? (*nodes_type)(local_id) : _nt_normal;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNodeGlobalId(UInt local_id) const {
  return nodes_global_ids ? (*nodes_global_ids)(local_id) : local_id;
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNodeLocalId(UInt global_id) const {
  AKANTU_DEBUG_ASSERT(nodes_global_ids != nullptr,
                      "The global ids are note set.");
  return nodes_global_ids->find(global_id);
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbGlobalNodes() const {
  return nodes_global_ids ? nb_global_nodes : nodes->size();
}

/* -------------------------------------------------------------------------- */
inline UInt Mesh::getNbNodesPerElementList(const Array<Element> & elements) {
  UInt nb_nodes_per_element = 0;
  UInt nb_nodes = 0;
  ElementType current_element_type = _not_defined;

  for(const auto & el : elements) {
    if (el.type != current_element_type) {
      current_element_type = el.type;
      nb_nodes_per_element = Mesh::getNbNodesPerElement(current_element_type);
    }

    nb_nodes += nb_nodes_per_element;
  }

  return nb_nodes;
}

/* -------------------------------------------------------------------------- */
inline Mesh & Mesh::getMeshFacets() {
  if (!this->mesh_facets)
    AKANTU_SILENT_EXCEPTION(
        "No facet mesh is defined yet! check the buildFacets functions");

  return *this->mesh_facets;
}

/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getMeshFacets() const {
  if (!this->mesh_facets)
    AKANTU_SILENT_EXCEPTION(
        "No facet mesh is defined yet! check the buildFacets functions");

  return *this->mesh_facets;
}
/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getMeshParent() const {
  if (!this->mesh_parent)
    AKANTU_SILENT_EXCEPTION(
        "No parent mesh is defined! This is only valid in a mesh_facets");

  return *this->mesh_parent;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_MESH_INLINE_IMPL_CC__ */
