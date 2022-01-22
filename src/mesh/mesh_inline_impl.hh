/**
 * @file   mesh_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Thu Jul 15 2010
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Implementation of the inline functions of the mesh class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_iterators.hh"
#include "element_class.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

// #ifndef __AKANTU_MESH_INLINE_IMPL_CC__
// #define __AKANTU_MESH_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getNbFacetsPerElement(ElementType type) -> Int {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getNbFacetsPerElement();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getNbFacetsPerElement(ElementType type, Idx t)
    -> Int {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getNbFacetsPerElement(t);
      },
      type);
}

/* -------------------------------------------------------------------------- */
template <typename... pack>
auto Mesh::elementTypes(pack &&... _pack) const -> ElementTypesIteratorHelper {
  return connectivities.elementTypes(_pack...);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getConnectivity(const Element & element) const {
  return connectivities.get(element);
}

/* -------------------------------------------------------------------------- */
inline RemovedNodesEvent::RemovedNodesEvent(const Mesh & mesh,
                                            const std::string & origin)
    : MeshEvent<Idx>(origin),
      new_numbering(mesh.getNbNodes(), 1, "new_numbering") {}

/* -------------------------------------------------------------------------- */
inline RemovedElementsEvent::RemovedElementsEvent(const Mesh & mesh,
                                                  const ID & new_numbering_id,
                                                  const std::string & origin)
    : MeshEvent<Element>(origin),
      new_numbering(new_numbering_id, mesh.getID()) {}

/* -------------------------------------------------------------------------- */
template <>
inline void Mesh::sendEvent<NewElementsEvent>(NewElementsEvent & event) {
  this->fillNodesToElements();
  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <> inline void Mesh::sendEvent<NewNodesEvent>(NewNodesEvent & event) {
  this->computeBoundingBox();
  this->nodes_flags->resize(this->nodes->size(), NodeFlag::_normal);
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
  if (nodes_global_ids and not is_mesh_facets) {
    this->removeNodesFromArray(*nodes_global_ids, new_numbering);
  }
  if (not is_mesh_facets) {
    this->removeNodesFromArray(*nodes_flags, new_numbering);
  }

  if (not nodes_to_elements.empty()) {
    std::vector<std::unique_ptr<std::set<Element>>> tmp(
        nodes_to_elements.size());
    auto it = nodes_to_elements.begin();

    Int new_nb_nodes = 0;
    for (auto new_i : new_numbering) {
      if (new_i != Int(-1)) {
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
                                       const Array<Idx> & new_numbering) {
  Array<T> tmp(vect.size(), vect.getNbComponent());
  auto nb_component = vect.getNbComponent();
  auto new_nb_nodes = 0;
  for (Int i = 0; i < new_numbering.size(); ++i) {
    auto new_i = new_numbering(i);
    if (new_i != Int(-1)) {
      T * to_copy = vect.data() + i * nb_component;
      std::uninitialized_copy(to_copy, to_copy + nb_component,
                              tmp.data() + new_i * nb_component);
      ++new_nb_nodes;
    }
  }

  tmp.resize(new_nb_nodes);
  vect.copy(tmp);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNodesGlobalIdsPointer() -> Array<Idx> & {
  AKANTU_DEBUG_IN();
  if (not nodes_global_ids) {
    nodes_global_ids = std::make_shared<Array<Idx>>(
        nodes->size(), 1, getID() + ":nodes_global_ids");

    for (auto && global_ids : enumerate(*nodes_global_ids)) {
      std::get<1>(global_ids) = std::get<0>(global_ids);
    }
  }

  AKANTU_DEBUG_OUT();
  return *nodes_global_ids;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto)
Mesh::getDataPointer(const ID & data_name, ElementType el_type,
                     GhostType ghost_type, Int nb_component,
                     bool size_to_nb_element, bool resize_with_parent) {
  Array<T> & tmp = this->getElementalDataArrayAlloc<T>(
      data_name, el_type, ghost_type, nb_component);

  if (size_to_nb_element) {
    if (resize_with_parent) {
      tmp.resize(mesh_parent->getNbElement(el_type, ghost_type));
    } else {
      tmp.resize(this->getNbElement(el_type, ghost_type));
    }
  }

  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto)
Mesh::getDataPointer(const ID & data_name, ElementType el_type,
                     GhostType ghost_type, Int nb_component,
                     bool size_to_nb_element, bool resize_with_parent,
                     const T & defaul_) {
  Array<T> & tmp = this->getElementalDataArrayAlloc<T>(
      data_name, el_type, ghost_type, nb_component);

  if (size_to_nb_element) {
    if (resize_with_parent) {
      tmp.resize(mesh_parent->getNbElement(el_type, ghost_type), defaul_);
    } else {
      tmp.resize(this->getNbElement(el_type, ghost_type), defaul_);
    }
  }

  return tmp;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto) Mesh::getData(const ID & data_name, ElementType el_type,
                                    GhostType ghost_type) const {
  return this->getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto) Mesh::getData(const ID & data_name, ElementType el_type,
                                    GhostType ghost_type) {
  return this->getElementalDataArray<T>(data_name, el_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto) Mesh::getData(const ID & data_name) const {
  return this->getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline decltype(auto) Mesh::getData(const ID & data_name) {
  return this->getElementalData<T>(data_name);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getConnectivityPointer(ElementType type, GhostType ghost_type)
    -> Array<Idx> & {
  if (connectivities.exists(type, ghost_type)) {
    return connectivities(type, ghost_type);
  }

  if (ghost_type != _not_ghost) {
    ghosts_counters.alloc(0, 1, type, ghost_type, 1);
  }

  AKANTU_DEBUG_INFO("The connectivity vector for the type " << type
                                                            << " created");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  return connectivities.alloc(0, nb_nodes_per_element, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getElementToSubelementPointer(ElementType type, GhostType ghost_type) {
  return getDataPointer<std::vector<Element>>("element_to_subelement", type,
                                              ghost_type, 1, true);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getSubelementToElementPointer(ElementType type, GhostType ghost_type) {
  auto & array = getDataPointer<Element>(
      "subelement_to_element", type, ghost_type, getNbFacetsPerElement(type),
      false, is_mesh_facets, ElementNull);
  return array;
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getElementToSubelement() const {
  return getData<std::vector<Element>>("element_to_subelement");
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getElementToSubelementNC() {
  return getData<std::vector<Element>>("element_to_subelement");
}

/* -------------------------------------------------------------------------- */
inline const auto & Mesh::getElementToSubelement(ElementType type,
                                                 GhostType ghost_type) const {
  return getData<std::vector<Element>>("element_to_subelement", type,
                                       ghost_type);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getElementToSubelementNC(ElementType type,
                                             GhostType ghost_type) {
  return getData<std::vector<Element>>("element_to_subelement", type,
                                       ghost_type);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getElementToSubelement(const Element & element) const {
  return getData<std::vector<Element>>("element_to_subelement")(element, 0);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getElementToSubelementNC(const Element & element) {
  return getData<std::vector<Element>>("element_to_subelement")(element, 0);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getSubelementToElement() const {
  return getData<Element>("subelement_to_element");
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getSubelementToElementNC() {
  return getData<Element>("subelement_to_element");
}

/* -------------------------------------------------------------------------- */
inline const auto & Mesh::getSubelementToElement(ElementType type,
                                                 GhostType ghost_type) const {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline auto & Mesh::getSubelementToElementNC(ElementType type,
                                             GhostType ghost_type) {
  return getData<Element>("subelement_to_element", type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getSubelementToElement(const Element & element) const {
  return this->getSubelementToElement().get(element);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getSubelementToElementNC(const Element & element) const {
  return this->getSubelementToElement().get(element);
}

/* -------------------------------------------------------------------------- */
template <class D, std::enable_if_t<aka::is_vector<D>::value> *>
inline void Mesh::getBarycenter(const Element & element,
                                Eigen::MatrixBase<D> & barycenter) const {
  const auto && conn = getConnectivity(element);
  Matrix<Real> local_coord(spatial_dimension, conn.size());
  auto node_begin = make_view(*nodes, spatial_dimension).begin();

  for (auto && data : enumerate(conn)) {
    local_coord(std::get<0>(data)) = node_begin[std::get<1>(data)];
  }

  Math::barycenter(local_coord, barycenter);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getKind(ElementType type) -> ElementKind {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getKind();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Element::kind() const -> ElementKind {
  return Mesh::getKind(type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getP1ElementType(ElementType type) -> ElementType {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getP1ElementType();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getSpatialDimension(ElementType type) -> Int {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getSpatialDimension();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getNaturalSpaceDimension(ElementType type) -> Int {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getNaturalSpaceDimension();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getNbFacetTypes(ElementType type, Idx /*t*/)
    -> Int {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getNbFacetTypes();
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline constexpr auto Mesh::getFacetType(ElementType type, Idx t)
    -> ElementType {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getFacetType(t);
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getAllFacetTypes(ElementType type) {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        auto && map = ElementClass<type>::getFacetTypes();
        return Eigen::Map<const Eigen::Matrix<ElementType, Eigen::Dynamic, 1>>(
            map.data(), map.rows(), map.cols());
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getFacetLocalConnectivity(ElementType type, Idx t) {
  return tuple_dispatch<AllElementTypes>(
      [&](auto && enum_type) {
        constexpr ElementType type = std::decay_t<decltype(enum_type)>::value;
        return ElementClass<type>::getFacetLocalConnectivityPerElement(t);
      },
      type);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getFacetConnectivity(const Element & element, Idx t) const
    -> Matrix<Idx> {
  auto local_facets = getFacetLocalConnectivity(element.type, t);
  Matrix<Idx> facets(local_facets.rows(), local_facets.cols());

  const auto & conn = connectivities(element.type, element.ghost_type);

  for (Int f = 0; f < facets.rows(); ++f) {
    for (Int n = 0; n < facets.cols(); ++n) {
      facets(f, n) = conn(element.element, local_facets(f, n));
    }
  }

  return facets;
}

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getConnectivityNC(const Element & element) {
  return connectivities.get(element);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Mesh::extractNodalValuesFromElement(
    const Array<T> & nodal_values, T * local_coord, Int * connectivity,
    Int n_nodes, Int nb_degree_of_freedom) const {
  for (Int n = 0; n < n_nodes; ++n) {
    memcpy(local_coord + n * nb_degree_of_freedom,
           nodal_values.data() + connectivity[n] * nb_degree_of_freedom,
           nb_degree_of_freedom * sizeof(T));
  }
}

/* -------------------------------------------------------------------------- */
inline void Mesh::addConnectivityType(ElementType type, GhostType ghost_type) {
  getConnectivityPointer(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isPureGhostNode(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_shared_mask) == NodeFlag::_pure_ghost;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isLocalOrMasterNode(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_local_master_mask) == NodeFlag::_normal;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isLocalNode(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_shared_mask) == NodeFlag::_normal;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isMasterNode(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_shared_mask) == NodeFlag::_master;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isSlaveNode(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_shared_mask) == NodeFlag::_slave;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isPeriodicSlave(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_periodic_mask) ==
         NodeFlag::_periodic_slave;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::isPeriodicMaster(Idx n) const -> bool {
  return ((*nodes_flags)(n)&NodeFlag::_periodic_mask) ==
         NodeFlag::_periodic_master;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNodeFlag(Idx local_id) const -> NodeFlag {
  return (*nodes_flags)(local_id);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNodePrank(Idx local_id) const {
  auto it = nodes_prank.find(local_id);
  return it == nodes_prank.end() ? -1 : it->second;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNodeGlobalId(Idx local_id) const {
  return nodes_global_ids ? (*nodes_global_ids)(local_id) : local_id;
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNodeLocalId(Idx global_id) const {
  if (nodes_global_ids == nullptr) {
    return global_id;
  }
  return nodes_global_ids->find(global_id);
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNbGlobalNodes() const {
  return nodes_global_ids ? nb_global_nodes : nodes->size();
}

/* -------------------------------------------------------------------------- */
inline auto Mesh::getNbNodesPerElementList(const Array<Element> & elements)
    -> Int {
  Int nb_nodes_per_element = 0;
  Int nb_nodes = 0;
  ElementType current_element_type = _not_defined;

  for (const auto & el : elements) {
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
  if (this->mesh_facets == nullptr) {
    AKANTU_SILENT_EXCEPTION(
        "No facet mesh is defined yet! check the buildFacets functions");
  }

  return *this->mesh_facets;
}

/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getMeshFacets() const {
  if (this->mesh_facets == nullptr) {
    AKANTU_SILENT_EXCEPTION(
        "No facet mesh is defined yet! check the buildFacets functions");
  }

  return *this->mesh_facets;
}
/* -------------------------------------------------------------------------- */
inline const Mesh & Mesh::getMeshParent() const {
  if (this->mesh_parent == nullptr) {
    AKANTU_SILENT_EXCEPTION(
        "No parent mesh is defined! This is only valid in a mesh_facets");
  }

  return *this->mesh_parent;
}

/* -------------------------------------------------------------------------- */
void Mesh::addPeriodicSlave(Idx slave, Idx master) {
  if (master == slave) {
    return;
  }

  // if pair already registered
  auto master_slaves = periodic_master_slave.equal_range(master);
  auto slave_it =
      std::find_if(master_slaves.first, master_slaves.second,
                   [&](auto & pair) { return pair.second == slave; });
  if (slave_it == master_slaves.second) {
    // no duplicates
    periodic_master_slave.insert(std::make_pair(master, slave));
    AKANTU_DEBUG_INFO("adding periodic slave, slave gid:"
                      << getNodeGlobalId(slave) << " [lid: " << slave << "]"
                      << ", master gid:" << getNodeGlobalId(master)
                      << " [lid: " << master << "]");
    // std::cout << "adding periodic slave, slave gid:" <<
    // getNodeGlobalId(slave)
    //           << " [lid: " << slave << "]"
    //           << ", master gid:" << getNodeGlobalId(master)
    //           << " [lid: " << master << "]" << std::endl;
  }

  periodic_slave_master[slave] = master;

  auto set_flag = [&](auto node, auto flag) {
    (*nodes_flags)[node] &= ~NodeFlag::_periodic_mask; // clean periodic flags
    (*nodes_flags)[node] |= flag;
  };

  set_flag(slave, NodeFlag::_periodic_slave);
  set_flag(master, NodeFlag::_periodic_master);
}

/* -------------------------------------------------------------------------- */
auto Mesh::getPeriodicMaster(Idx slave) const -> Idx {
  return periodic_slave_master.at(slave);
}

/* -------------------------------------------------------------------------- */
class Mesh::PeriodicSlaves {
  using internal_iterator = std::unordered_multimap<Idx, Idx>::const_iterator;
  std::pair<internal_iterator, internal_iterator> pair;

public:
  PeriodicSlaves(const Mesh & mesh, Idx master)
      : pair(mesh.getPeriodicMasterSlaves().equal_range(master)) {}

  PeriodicSlaves(const PeriodicSlaves & other) = default;
  PeriodicSlaves(PeriodicSlaves && other) noexcept = default;
  auto operator=(const PeriodicSlaves & other) -> PeriodicSlaves & = default;

  class const_iterator {
    internal_iterator it;

  public:
    const_iterator(internal_iterator it) : it(it) {}

    const_iterator operator++() {
      ++it;
      return *this;
    }
    bool operator!=(const const_iterator & other) { return other.it != it; }
    bool operator==(const const_iterator & other) { return other.it == it; }
    auto operator*() { return it->second; }
  };

  auto begin() const { return const_iterator(pair.first); }
  auto end() const { return const_iterator(pair.second); }
};

/* -------------------------------------------------------------------------- */
inline decltype(auto) Mesh::getPeriodicSlaves(Idx master) const {
  return PeriodicSlaves(*this, master);
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
Mesh::getConnectivityWithPeriodicity(const Element & element) const {
  Vector<Idx> conn = getConnectivity(element);
  if (not isPeriodic()) {
    return conn;
  }

  for (auto && node : conn) {
    if (isPeriodicSlave(node)) {
      node = getPeriodicMaster(node);
    }
  }

  return conn;
}

} // namespace akantu

//#endif /* __AKANTU_MESH_INLINE_IMPL_CC__ */
