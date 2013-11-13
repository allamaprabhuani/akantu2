/**
 * @file   node_group_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Jun  7 12:57:24 2013
 *
 * @brief  Node group inline function definitions
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::begin() const {
  return node_group.begin();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::end() const {
  return node_group.end();
}

/* -------------------------------------------------------------------------- */
inline NodeGroup::const_node_iterator NodeGroup::add(UInt node) {
  const_node_iterator it = std::find(begin(), end(), node);
  if(it == node_group.end()) {
    node_group.push_back(node);
    return (node_group.end() - 1);
  }

  return it;
}

/* -------------------------------------------------------------------------- */
inline UInt NodeGroup::getSize() const {
  return node_group.getSize();
}


__END_AKANTU__
