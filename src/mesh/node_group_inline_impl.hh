/**
 * @file   node_group_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Dec 09 2020
 *
 * @brief  Node group inline function definitions
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
#include "node_group.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::begin() const { return node_group.begin(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::end() const { return node_group.end(); }

/* -------------------------------------------------------------------------- */
inline auto NodeGroup::add(Idx node, bool check_for_duplicate) {
  const_node_iterator it;
  if (check_for_duplicate) {
    it = std::find(begin(), end(), node);
    if (it != node_group.end()) {
      return it;
    }
  }

  node_group.push_back(node);
  it = (node_group.end() - 1);
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

} // namespace akantu
