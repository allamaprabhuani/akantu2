/**
 * @file   mesh_utils.cc
 * @author Guillaume ANCIAUX <anciaux@lsmscluster1.epfl.ch>
 * @date   Wed Aug 18 14:21:00 2010
 *
 * @brief  All mesh utils necessary for various tasks
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#include "mesh_utils.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2Elements(const Mesh & mesh,
				   Vector<UInt> & node_offset,
				   Vector<UInt> & node_to_elem) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  ElementType type_p1;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != mesh.getSpatialDimension()) continue;

    nb_nodes_per_element[nb_good_types]    = Mesh::getNbNodesPerElement(type);
    type_p1 = Mesh::getP1ElementType(type);
    nb_nodes_per_element_p1[nb_good_types] = Mesh::getNbNodesPerElement(type_p1);

    conn_val[nb_good_types] = mesh.getConnectivity(type).values;
    nb_element[nb_good_types] = mesh.getConnectivity(type).getSize();
    nb_good_types++;
  }

  AKANTU_DEBUG_ASSERT(nb_good_types  != 0,
		      "Some elements must be found in right dimension to compute facets!");

  /// array for the node-element list
  node_offset.resize(nb_nodes + 1);
  UInt * node_offset_val = node_offset.values;

  /// count number of occurrence of each node
  for (UInt t = 0; t < nb_good_types; ++t) {
    memset(node_offset_val, 0, (nb_nodes + 1)*sizeof(UInt));
    for (UInt el = 0; el < nb_element[t]; ++el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n) {
	node_offset_val[conn_val[t][el_offset + n]]++;
      }
    }
  }

  /// convert the occurrence array in a csr one
  for (UInt i = 1; i < nb_nodes; ++i) node_offset_val[i] += node_offset_val[i-1];
  for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  node_offset_val[0] = 0;

  /// rearrange element to get the node-element list
  node_to_elem.resize(node_offset_val[nb_nodes]);
  UInt * node_to_elem_val = node_to_elem.values;

  for (UInt t = 0, linearized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linearized_el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element_p1[t]; ++n)
	node_to_elem_val[node_offset_val[conn_val[t][el_offset + n]]++] = linearized_el;
    }

  for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  node_offset_val[0] = 0;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacets(Mesh & mesh, bool boundary_flag, bool internal_flag){
  AKANTU_DEBUG_IN();

  Vector<UInt> node_offset;
  Vector<UInt> node_to_elem;

  buildNode2Elements(mesh,node_offset,node_to_elem);

  //  std::cout << "node offset " << std::endl << node_offset << std::endl;
  // std::cout << "node to elem " << std::endl << node_to_elem << std::endl;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_facet[nb_types];

  UInt nb_facets[nb_types];
  UInt ** node_in_facet[nb_types];
  Vector<UInt> * connectivity_facets[nb_types];
  Vector<UInt> * connectivity_internal_facets[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  ElementType facet_type;

  Mesh * internal_facets_mesh = NULL;
  UInt spatial_dimension = mesh.getSpatialDimension();
  if (internal_flag) internal_facets_mesh = mesh.getInternalFacetsMeshPointer();

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(mesh.getSpatialDimension(type) != spatial_dimension) continue;

    nb_nodes_per_element[nb_good_types]    = mesh.getNbNodesPerElement(type);

    facet_type = mesh.getFacetElementType(type);
    nb_facets[nb_good_types] = mesh.getNbFacetsPerElementType(type);
    node_in_facet[nb_good_types] = mesh.getFacetLocalConnectivityPerElementType(type);

    nb_nodes_per_facet[nb_good_types]    = mesh.getNbNodesPerElement(facet_type);

    if (boundary_flag){
      // getting connectivity of boundary facets
      connectivity_facets[nb_good_types] = mesh.getConnectivityPointer(facet_type);
      connectivity_facets[nb_good_types]->resize(0);
    }
    if (internal_flag){
      // getting connectivity of internal facets
      connectivity_internal_facets[nb_good_types] =
	internal_facets_mesh->getConnectivityPointer(facet_type);
      connectivity_internal_facets[nb_good_types]->resize(0);
    }

    conn_val[nb_good_types] = mesh.getConnectivity(type).values;
    nb_element[nb_good_types] = mesh.getConnectivity(type).getSize();
    nb_good_types++;
  }

  Vector<UInt> counter;
  /// count number of occurrence of each node
  for (UInt t = 0,linearized_el = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el, ++linearized_el) {

      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt f = 0; f < nb_facets[t]; ++f) {
	//build the nodes involved in facet 'f'
	UInt facet_nodes[nb_nodes_per_facet[t]];
	for (UInt n = 0; n < nb_nodes_per_facet[t]; ++n) {
	  UInt node_facet = node_in_facet[t][f][n];
	  facet_nodes[n] = conn_val[t][el_offset + node_facet];
	}
	//our reference is the first node
	UInt * first_node_elements = node_to_elem.values+node_offset.values[facet_nodes[0]];
	UInt first_node_nelements =
	  node_offset.values[facet_nodes[0]+1]-
	  node_offset.values[facet_nodes[0]];
	counter.resize(first_node_nelements);
	memset(counter.values,0,sizeof(UInt)*first_node_nelements);
	//loop over the other nodes to search intersecting elements
	for (UInt n = 1; n < nb_nodes_per_facet[t]; ++n) {
	  UInt * node_elements = node_to_elem.values+node_offset.values[facet_nodes[n]];
	  UInt node_nelements =
	    node_offset.values[facet_nodes[n]+1]-
	    node_offset.values[facet_nodes[n]];
	  for (UInt el1 = 0; el1 < first_node_nelements; ++el1) {
	    for (UInt el2 = 0; el2 < node_nelements; ++el2) {
	      if (first_node_elements[el1] == node_elements[el2]) {
		++counter.values[el1];
		break;
	      }
	    }
	    if (counter.values[el1] == nb_nodes_per_facet[t]) break;
	  }
	}
	bool connected_facet = false;
	//the connected elements are those for which counter is n_facet
	//	UInt connected_element = -1;
	for (UInt el1 = 0; el1 < counter.getSize(); el1++) {
	  UInt el_index = node_to_elem.values[node_offset.values[facet_nodes[0]]+el1];
	  if (counter.values[el1] == nb_nodes_per_facet[t]-1 && el_index > linearized_el){
	    //	    connected_element = el_index;
	    AKANTU_DEBUG(dblDump,"connecting elements " << linearized_el << " and " << el_index);   
	    if (internal_flag)
	      connectivity_internal_facets[t]->push_back(facet_nodes);
	  }
	  if (counter.values[el1] == nb_nodes_per_facet[t]-1 && el_index != linearized_el)
	    connected_facet = true;
	}
	if (!connected_facet) {
	  AKANTU_DEBUG(dblDump,"facet " << f << " in element " << linearized_el << " is a boundary");
	  if (boundary_flag)
	    connectivity_facets[t]->push_back(facet_nodes);
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNormals(Mesh & mesh,UInt spatial_dimension){
  AKANTU_DEBUG_IN();
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;


  UInt nb_types = type_list.size();
  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];
  UInt nb_good_types = 0;

  Vector<UInt> * connectivity[nb_types];
  Vector<Real> * normals[nb_types];

  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    ElementType type_p1 = mesh.getP1ElementType(type);
    if(mesh.getSpatialDimension(type) != spatial_dimension) continue;

    nb_nodes_per_element[nb_good_types]    = mesh.getNbNodesPerElement(type);
    nb_nodes_per_element_p1[nb_good_types] = mesh.getNbNodesPerElement(type_p1);

    // getting connectivity
    connectivity[nb_good_types] = mesh.getConnectivityPointer(type);
    if (!connectivity[nb_good_types])
      AKANTU_DEBUG_ERROR("connectivity is not allocatted : this should probably not have happend");

    //getting array of normals
    normals[nb_good_types] = mesh.getNormalsPointer(type);
    if(normals[nb_good_types])
      normals[nb_good_types]->resize(0);
    else
      normals[nb_good_types] = &mesh.createNormals(type);


    nb_good_types++;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::renumberMeshNodes(Mesh & mesh,
				  UInt * local_connectivities,
				  UInt nb_local_element,
				  UInt nb_ghost_element,
				  ElementType type,
				  Vector<UInt> & nodes_numbers) {
  AKANTU_DEBUG_IN();

  nodes_numbers.resize(0);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  UInt * conn = local_connectivities;

  /// renumber the nodes
  for (UInt el = 0; el < nb_local_element + nb_ghost_element; ++el) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Int nn = nodes_numbers.find(*conn);
      if(nn == -1) {
	nodes_numbers.push_back(*conn);
	*conn = nodes_numbers.getSize() - 1;
      } else {
	*conn = nn;
      }
      *conn++;
    }
  }

  /// copy the renumbered connectivity to the right place
  Vector<UInt> * local_conn = mesh.getConnectivityPointer(type);
  local_conn->resize(nb_local_element);
  memcpy(local_conn->values,
	 local_connectivities,
	 nb_local_element * nb_nodes_per_element * sizeof(UInt));

  Vector<UInt> * ghost_conn = mesh.getGhostConnectivityPointer(type);
  ghost_conn->resize(nb_ghost_element);
  memcpy(ghost_conn->values,
	 local_connectivities + nb_local_element * nb_nodes_per_element,
	 nb_ghost_element * nb_nodes_per_element * sizeof(UInt));

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
