/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "aka_iterators.hh"
#include "element_group.hh"
#include "element_type_map.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

//#ifndef AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_
//#define AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void ElementGroup::add(const Element & el, bool add_nodes,
                              bool check_for_duplicate) {
  this->add(el.type, el.element, el.ghost_type, add_nodes, check_for_duplicate);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::add(ElementType type, Idx element,
                              GhostType ghost_type, bool add_nodes,
                              bool check_for_duplicate) {
  addElement(type, element, ghost_type);

  if (add_nodes) {
    auto it = mesh.getConnectivity(type, ghost_type)
                  .begin(mesh.getNbNodesPerElement(type)) +
              element;
    auto && conn = *it;
    for (Idx i = 0; i < conn.size(); ++i) {
      addNode(conn[i], check_for_duplicate);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addNode(Idx node_id, bool check_for_duplicate) {
  node_group.add(node_id, check_for_duplicate);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::removeNode(Idx node_id) {
  node_group.remove(node_id);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addElement(ElementType elem_type, Idx elem_id,
                                     GhostType ghost_type) {
  if (!(elements.exists(elem_type, ghost_type))) {
    elements.alloc(0, 1, elem_type, ghost_type);
  }

  elements(elem_type, ghost_type).push_back(elem_id);
  this->dimension = Int(
      std::max(Int(this->dimension), Int(mesh.getSpatialDimension(elem_type))));
}

/* -------------------------------------------------------------------------- */
inline Int ElementGroup::getNbNodes() const { return node_group.size(); }

/* -------------------------------------------------------------------------- */
inline auto ElementGroup::begin(ElementType type, GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type).begin();
  }
  return empty_elements.begin();
}

/* -------------------------------------------------------------------------- */
inline auto ElementGroup::end(ElementType type, GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type).end();
  }
  return empty_elements.end();
}

/* -------------------------------------------------------------------------- */
inline const Array<Idx> &
ElementGroup::getElements(ElementType type, GhostType ghost_type) const {
  if (elements.exists(type, ghost_type)) {
    return elements(type, ghost_type);
  }
  return empty_elements;
}

/* -------------------------------------------------------------------------- */
inline decltype(auto)
ElementGroup::getElementsIterable(ElementType type,
                                  GhostType ghost_type) const {
  return make_transform_adaptor(this->elements(type, ghost_type),
                                [type, ghost_type](auto && el) {
                                  return Element{type, el, ghost_type};
                                });
}

} // namespace akantu

//#endif /* AKANTU_ELEMENT_GROUP_INLINE_IMPL_HH_ */
