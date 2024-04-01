/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "node_group.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::begin() const { return node_group.begin(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::end() const { return node_group.end(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::cbegin() const { return node_group.cbegin(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::cend() const { return node_group.cend(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::add(Idx node, bool check_for_duplicate) {
  if (check_for_duplicate) {
    auto it = std::find(cbegin(), cend(), node);
    if (it != node_group.cend()) {
      return it;
    }
  }

  node_group.push_back(node);
  auto it = (node_group.cend() - 1);
  return it;
}

/* -------------------------------------------------------------------------- */
inline void NodeGroup::remove(Idx node) {
  auto it = this->node_group.begin();
  auto end = this->node_group.end();
  AKANTU_DEBUG_ASSERT(it != end, "The node group is empty!!");
  for (; it != node_group.end(); ++it) {
    if (*it == node) {
      it = node_group.erase(it);
    }
  }
  AKANTU_DEBUG_ASSERT(it != end, "The node was not found!");
}

/* -------------------------------------------------------------------------- */
inline bool NodeGroup::empty() const { return node_group.empty(); }

/* -------------------------------------------------------------------------- */
inline Int NodeGroup::size() const { return node_group.size(); }

/* -------------------------------------------------------------------------- */
struct FilterFunctor;

/* -------------------------------------------------------------------------- */
template <typename T> void NodeGroup::applyNodeFilter(T & filter) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(T::type == FilterFunctor::_node_filter_functor,
                      "NodeFilter can only apply node filter functor");

  auto it = this->node_group.begin();

  for (; it != node_group.end(); ++it) {
    /// filter == true -> keep node
    if (!filter(*it)) {
      it = node_group.erase(it);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::getNbGlobalNodes() const -> Idx {
  return nb_global_nodes == -1 ? node_group.size() : nb_global_nodes;
}

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::getNbNodes() const -> Idx { return node_group.size(); }

} // namespace akantu
