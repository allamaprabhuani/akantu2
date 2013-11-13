/**
 * @file   element_group_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information relevent to the notion of domain boundary and surfaces.
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
inline void ElementGroup::add(const Element & el) {
  addElement(el.type, el.element, el.ghost_type);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addNode(UInt node_id) {
  node_group.add(node_id);
}

/* -------------------------------------------------------------------------- */
inline void ElementGroup::addElement(const ElementType & elem_type,
				    UInt elem_id,
				    const GhostType & ghost_type) {
  if(!(elements.exists(elem_type, ghost_type))) {
    elements.alloc(0, 1, elem_type, ghost_type);
  }
  elements(elem_type, ghost_type).push_back(elem_id);
}

/* -------------------------------------------------------------------------- */
inline UInt ElementGroup::getNbNodes() const {
  return node_group.getSize();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::const_node_iterator ElementGroup::node_begin() const {
  return node_group.begin();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::const_node_iterator ElementGroup::node_end() const {
  return node_group.end();
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::type_iterator ElementGroup::firstType(UInt dim,
                                                           const GhostType & ghost_type,
                                                           const ElementKind & kind) const {
  return elements.firstType(dim, ghost_type, kind);
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::type_iterator  ElementGroup::lastType(UInt dim,
                                                           const GhostType & ghost_type,
                                                           const ElementKind & kind) const {
  return elements.lastType(dim, ghost_type, kind);
}

/* -------------------------------------------------------------------------- */
inline ElementGroup::const_element_iterator ElementGroup::element_begin(const ElementType & type,
                                                                        const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(elements.exists(type, ghost_type), "No element of type " << type << ":" << ghost_type << " in the group " << name);

  return elements(type, ghost_type).begin();
}
 
/* -------------------------------------------------------------------------- */
inline ElementGroup::const_element_iterator ElementGroup::element_end(const ElementType & type,
                                                                      const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(elements.exists(type, ghost_type), "No element of type " << type << ":" << ghost_type << " in the group " << name);

  return elements(type, ghost_type).end();
}
