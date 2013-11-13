/**
 * @file   group_manager.cc
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author David Kammer <david.kammer@epfl.ch>
 *
 * @date   Wed Mar 06 09:30:00 2013
 *
 * @brief  Stores information about ElementGroup and NodeGroup
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
#include "group_manager.hh"
#include "mesh.hh"
#include "aka_csr.hh"
#include "mesh_utils.hh"

#include "element_group.hh"
#include "node_group.hh"
/* -------------------------------------------------------------------------- */
#include <sstream>
#include <algorithm>
#include <iterator>
#include <list>
#include <queue>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
GroupManager::GroupManager(const Mesh & mesh,
                           const ID & id,
                           const MemoryID & mem_id) : id(id),
						      memory_id(memory_id),
                                                      mesh(mesh) {

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
GroupManager::~GroupManager() {
  ElementGroups::iterator eit = element_groups.begin();
  ElementGroups::iterator eend = element_groups.end();
  for(; eit != eend ; ++eit) delete (eit->second);

  NodeGroups::iterator nit = node_groups.begin();
  NodeGroups::iterator nend = node_groups.end();
  for(; nit != nend ; ++nit) delete (nit->second);
}

/* -------------------------------------------------------------------------- */
NodeGroup & GroupManager::createNodeGroup(const std::string & group_name) {
  AKANTU_DEBUG_IN();

  if(node_groups.find(group_name + "_nodes") != node_groups.end()) {
    AKANTU_EXCEPTION("Trying to create a node group that already exists:" << group_name << "_nodes");
  }

  NodeGroup * node_group = new NodeGroup(group_name, id + ":" + group_name + "_node_group",
                                         memory_id);

  node_groups[group_name] = node_group;

  AKANTU_DEBUG_OUT();

  return *node_group;
}


/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name, UInt dimension) {
  AKANTU_DEBUG_IN();
  if(node_groups.find(group_name + "_nodes") != node_groups.end()) {
    AKANTU_EXCEPTION("Trying to create a node group that already exists:" << group_name << "_nodes");
  }

  if(element_groups.find(group_name) != element_groups.end()) {
    AKANTU_EXCEPTION("Trying to create a element group that already exists:" << group_name);
  }

  NodeGroup * node_group = new NodeGroup(group_name + "_nodes", id + ":" + group_name + "_node_group",
                                         memory_id);

  ElementGroup * element_group = new ElementGroup(group_name, mesh, *node_group,
                                                  dimension,
                                                  id + ":" + group_name + "_element_group",
                                                  memory_id);

  node_groups[group_name + "_nodes"] = node_group;
  element_groups[group_name] = element_group;

  AKANTU_DEBUG_OUT();

  return *element_group;
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
                                                UInt dimension,
                                                NodeGroup & node_group) {
  AKANTU_DEBUG_IN();

  if(element_groups.find(group_name) != element_groups.end()) {
    AKANTU_EXCEPTION("Trying to create a element group that already exists:" << group_name);
  }

  ElementGroup * element_group = new ElementGroup(group_name, mesh, node_group,
                                                  dimension,
                                                  id + ":" + group_name + "_element_group",
                                                  memory_id);

  element_groups[group_name] = element_group;

  AKANTU_DEBUG_OUT();

  return *element_group;
}

/* -------------------------------------------------------------------------- */
void GroupManager::createBoundaryGroupFromGeometry() {
  AKANTU_DEBUG_IN();

  CSR<Element> node_to_elem;

  /// Get list of surface elements
  UInt spatial_dimension = mesh.getSpatialDimension();

  //  buildNode2Elements(mesh, node_offset, node_to_elem, spatial_dimension-1);
  MeshUtils::buildNode2Elements(mesh, node_to_elem, spatial_dimension-1);

  std::list<Element> boundaries_elements;
  Mesh::type_iterator type_it = mesh.firstType(spatial_dimension - 1);
  Mesh::type_iterator type_end  = mesh.lastType(spatial_dimension - 1);
  Element el; el.ghost_type = _not_ghost;
  for (; type_it != type_end; ++type_it) {
    el.type = *type_it;
    for (UInt e = 0; e < mesh.getNbElement(*type_it, _not_ghost); ++e) {
      el.element = e;
      boundaries_elements.push_back(el);
    }
  }

  UInt nb_boundaries = 0;

  while(!boundaries_elements.empty()) {
    std::queue<Element> element_to_add;
    element_to_add.push(*(boundaries_elements.begin()));
    boundaries_elements.erase(boundaries_elements.begin());

    std::stringstream sstr; sstr << "boundary_" << nb_boundaries;
    ElementGroup & boundary = createElementGroup(sstr.str(), spatial_dimension - 1);
    ++nb_boundaries;

    while(!element_to_add.empty()) {
      Element el = element_to_add.front();
      boundary.add(el);
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(el.type);
      element_to_add.pop();
      Vector<UInt> connect =
	mesh.getConnectivity(el.type, el.ghost_type).begin(nb_nodes_per_element)[el.element];
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt node = connect[n];
	boundary.addNode(node);

	CSR<Element>::iterator it_n;
	for (it_n = node_to_elem.begin(node); it_n != node_to_elem.end(node); ++it_n) {
	  Element & check_el = *it_n;
	  std::list<Element>::iterator it_boun = std::find(boundaries_elements.begin(),
							   boundaries_elements.end(),
							   check_el);
	  if(it_boun != boundaries_elements.end()) {
	    boundaries_elements.erase(it_boun);
	    element_to_add.push(check_el);
	  }
	}
      }
    }
    boundary.getNodeGroup().removeDuplicate();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename T>
void GroupManager::createGroupsFromMeshData(const std::string & dataset_name) {
  std::set<std::string> group_names;
  const ByElementTypeArray<T> & datas = mesh.getData<T>(dataset_name);
  typedef typename ByElementTypeArray<T>::type_iterator type_iterator;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    type_iterator type_it = datas.firstType(_all_dimensions, *gt);
    type_iterator type_end  = datas.lastType(_all_dimensions, *gt);
    for (; type_it != type_end; ++type_it) {
      const Array<T> & dataset = datas(*type_it, *gt);
      UInt nb_element = mesh.getNbElement(*type_it, *gt);
      AKANTU_DEBUG_ASSERT(dataset.getSize() == nb_element,
			  "Not the same number of elements in the map from MeshData and in the mesh!");
      for(UInt e(0); e < nb_element; ++e) {
	std::stringstream sstr; sstr << dataset(e);
	group_names.insert(sstr.str());
      }
    }
  }

  /// @todo sharing the group names in parallel to avoid trouble with groups
  /// present only on a sub set of processors
  std::set<std::string>::iterator git  = group_names.begin();
  std::set<std::string>::iterator gend = group_names.end();
  for (;git != gend; ++git) createElementGroup(*git, _all_dimensions);

  Element el;
  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    el.ghost_type = *gt;

    type_iterator type_it = datas.firstType(_all_dimensions, *gt);
    type_iterator type_end  = datas.lastType(_all_dimensions, *gt);
    for (; type_it != type_end; ++type_it) {
      el.type = *type_it;

      const Array<T> & dataset = datas(*type_it, *gt);
      UInt nb_element = mesh.getNbElement(*type_it, *gt);
      AKANTU_DEBUG_ASSERT(dataset.getSize() == nb_element,
			  "Not the same number of elements in the map from MeshData and in the mesh!");

      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(el.type);

      Array<UInt>::const_iterator< Vector<UInt> > cit =
	mesh.getConnectivity(*type_it, *gt).begin(nb_nodes_per_element);

      for(UInt e(0); e < nb_element; ++e, ++cit) {
	el.element = e;
	std::stringstream sstr; sstr << dataset(e);
	ElementGroup & group = getElementGroup(sstr.str());
	group.add(el);

	const Vector<UInt> & connect = *cit;
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  UInt node = connect[n];
	  group.addNode(node);
	}

      }
    }
  }

  git  = group_names.begin();
  for (;git != gend; ++git) getElementGroup(*git).getNodeGroup().removeDuplicate();
 }

template void GroupManager::createGroupsFromMeshData<std::string>(const std::string & dataset_name);
template void GroupManager::createGroupsFromMeshData<UInt>(const std::string & dataset_name);

/* -------------------------------------------------------------------------- */
void GroupManager::createElementGroupFromNodeGroup(const std::string & name,
						   const std::string & node_group_name,
						   UInt dimension) {
  NodeGroup & node_group = getNodeGroup(node_group_name);
  ElementGroup & group = createElementGroup(name, dimension, node_group);

  CSR<Element> node_to_elem;
  MeshUtils::buildNode2Elements(mesh, node_to_elem, _all_dimensions);

  std::set<Element> seen;

  Array<UInt>::const_iterator<> itn  = node_group.begin();
  Array<UInt>::const_iterator<> endn = node_group.end();
  for (;itn != endn; ++itn) {
    CSR<Element>::iterator ite = node_to_elem.begin(*itn);
    CSR<Element>::iterator ende = node_to_elem.end(*itn);
    for (;ite != ende; ++ite) {
      const Element & elem = *ite;
      if(dimension != _all_dimensions && dimension != Mesh::getSpatialDimension(elem.type)) continue;
      if(seen.find(elem) != seen.end()) continue;

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(elem.type);
      Array<UInt>::const_iterator< Vector<UInt> > conn_it =
	mesh.getConnectivity(elem.type, elem.ghost_type).begin(nb_nodes_per_element);
      const Vector<UInt> & conn = conn_it[elem.element];

      UInt count = 0;
      for (UInt n = 0; n < conn.size(); ++n) {
	count += (node_group.getNodes().find(conn(n)) != -1 ? 1 : 0);
      }

      if(count == nb_nodes_per_element) group.add(elem);

      seen.insert(elem);
    }
  }
}

/* -------------------------------------------------------------------------- */
void GroupManager::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  stream << space << "GroupManager [" << std::endl;

  std::set<std::string> node_group_seen;
  for(const_element_group_iterator it(element_group_begin()); it != element_group_end(); ++it) {
    it->second->printself(stream, indent + 1);
    node_group_seen.insert(it->second->getNodeGroup().getName());
  }

  for(const_node_group_iterator it(node_group_begin()); it != node_group_end(); ++it) {
    if(node_group_seen.find(it->second->getName()) == node_group_seen.end())
      it->second->printself(stream, indent + 1);
  }

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
UInt GroupManager::getNbElementGroups(UInt dimension) const {
  if(dimension == _all_dimensions) return element_groups.size();

  ElementGroups::const_iterator it  = element_groups.begin();
  ElementGroups::const_iterator end = element_groups.end();
  UInt count = 0;
  for(;it != end; ++it) count += (it->second->getDimension() == dimension);
  return count;
}

__END_AKANTU__
