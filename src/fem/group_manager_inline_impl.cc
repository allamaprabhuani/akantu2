/**
 * @file   boundary_inline_impl.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline const ElementGroup & GroupManager::getElementGroup(const std::string & name) const {
  const_element_group_iterator it = element_group_find(name);
  if(it == element_group_end()) {
    AKANTU_EXCEPTION("There are no element groups named " << name
		     << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline ElementGroup & GroupManager::getElementGroup(const std::string & name) {
  element_group_iterator it = element_group_find(name);
  if(it == element_group_end()) {
    AKANTU_EXCEPTION("There are no element groups named " << name
		     << " associated to the group manager: " << id);
  }

  return *(it->second);
}


/* -------------------------------------------------------------------------- */
inline const NodeGroup & GroupManager::getNodeGroup(const std::string & name) const {
  const_node_group_iterator it = node_group_find(name);
  if(it == node_group_end()) {
    AKANTU_EXCEPTION("There are no node groups named " << name
		     << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
inline NodeGroup & GroupManager::getNodeGroup(const std::string & name) {
  node_group_iterator it = node_group_find(name);
  if(it == node_group_end()) {
    AKANTU_EXCEPTION("There are no node groups named " << name
		     << " associated to the group manager: " << id);
  }

  return *(it->second);
}

__END_AKANTU__

