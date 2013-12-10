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
						      memory_id(mem_id),
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
NodeGroup & GroupManager::createNodeGroup(const std::string & group_name,
					  bool replace_group) {
  AKANTU_DEBUG_IN();

  NodeGroups::iterator it = node_groups.find(group_name);

  if(it != node_groups.end()) {
    if (replace_group) {
      it->second->empty();
      AKANTU_DEBUG_OUT();
      return *(it->second);
    }
    else
      AKANTU_EXCEPTION("Trying to create a node group that already exists:" << group_name << "_nodes");
  }

  NodeGroup * node_group = new NodeGroup(group_name, id + ":" + group_name + "_node_group",
                                         memory_id);

  node_groups[group_name] = node_group;

  AKANTU_DEBUG_OUT();

  return *node_group;
}


/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
						UInt dimension,
						bool replace_group) {
  AKANTU_DEBUG_IN();

  NodeGroup & new_node_group = createNodeGroup(group_name + "_nodes", replace_group);

  ElementGroups::iterator it = element_groups.find(group_name);

  if(it != element_groups.end()) {
    if (replace_group) {
      it->second->empty();
      AKANTU_DEBUG_OUT();
      return *(it->second);
    }
    else
      AKANTU_EXCEPTION("Trying to create a element group that already exists:" << group_name);
  }

  ElementGroup * element_group = new ElementGroup(group_name, mesh, new_node_group,
                                                  dimension,
                                                  id + ":" + group_name + "_element_group",
                                                  memory_id);

  node_groups[group_name + "_nodes"] = &new_node_group;
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
UInt GroupManager::createBoundaryGroupFromGeometry() {
  UInt spatial_dimension = mesh.getSpatialDimension();
  return createClusters(spatial_dimension - 1, "boundary_");
}

/* -------------------------------------------------------------------------- */
//// \todo if needed element list construction can be optimized by
//// templating the filter class
UInt GroupManager::createClusters(UInt element_dimension,
				  std::string cluster_name_prefix,
				  const GroupManager::ClusteringFilter & filter) {
  AKANTU_DEBUG_IN();

  CSR<Element> node_to_elem;

  /// Get list of all elements to check
  MeshUtils::buildNode2Elements(mesh, node_to_elem, element_dimension);

  std::list<Element> unseen_elements;
  Element el;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;
    el.ghost_type = ghost_type;

    Mesh::type_iterator type_it  = mesh.firstType(element_dimension,
						  ghost_type, _ek_not_defined);
    Mesh::type_iterator type_end = mesh.lastType (element_dimension,
						  ghost_type, _ek_not_defined);

    for (; type_it != type_end; ++type_it) {
      el.type = *type_it;
      el.kind = Mesh::getKind(*type_it);
      for (UInt e = 0; e < mesh.getNbElement(*type_it, ghost_type); ++e) {
	el.element = e;
	if (filter(el))
	  unseen_elements.push_back(el);
      }
    }
  }

  UInt nb_cluster = 0;

  /// keep looping until all elements are seen
  while(!unseen_elements.empty()) {
    /// create a new cluster
    std::stringstream sstr;
    sstr << cluster_name_prefix << nb_cluster;
    ElementGroup & cluster = createElementGroup(sstr.str(),
						element_dimension,
						true);
    ++nb_cluster;

    /// initialize the queue of elements to check in the current cluster
    std::queue<Element> element_to_add;
    element_to_add.push(*(unseen_elements.begin()));
    unseen_elements.erase(unseen_elements.begin());

    /// keep looping until current cluster is complete (no more
    /// connected elements)
    while(!element_to_add.empty()) {

      /// take first element and erase it in the queue
      Element el = element_to_add.front();
      element_to_add.pop();

      /// add current element to the cluster
      cluster.add(el);

      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
      Vector<UInt> connect =
	mesh.getConnectivity(el.type, el.ghost_type).begin(nb_nodes_per_element)[el.element];
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {

	/// add element's nodes to the cluster
	UInt node = connect[n];
	cluster.addNode(node);

	/// loop over neighbors
	CSR<Element>::iterator it_n;
	for (it_n = node_to_elem.begin(node); it_n != node_to_elem.end(node); ++it_n) {
	  Element & check_el = *it_n;
	  std::list<Element>::iterator it_clus = std::find(unseen_elements.begin(),
							   unseen_elements.end(),
							   check_el);

	  /// if neighbor not seen yet, add it to check list and
	  /// remove it from unseen elements
	  if(it_clus != unseen_elements.end()) {
	    unseen_elements.erase(it_clus);
	    element_to_add.push(check_el);
	  }
	}
      }
    }
    /// remove duplicated nodes
    cluster.getNodeGroup().removeDuplicate();
  }

  AKANTU_DEBUG_OUT();
  return nb_cluster;
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
