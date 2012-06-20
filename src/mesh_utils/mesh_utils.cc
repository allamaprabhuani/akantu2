/**
 * @file   mesh_utils.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Wed Aug 18 14:21:00 2010
 *
 * @brief  All mesh utils necessary for various tasks
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

#include "mesh_utils.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2Elements(const Mesh & mesh,
				   CSR<UInt> & node_to_elem,
				   UInt spatial_dimension) {
  AKANTU_DEBUG_IN();
  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes = mesh.getNbNodes();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;

  UInt nb_nodes_per_element[nb_types];
  //  UInt nb_nodes_per_element_p1[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  //  ElementType type_p1;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    nb_nodes_per_element[nb_good_types]    = Mesh::getNbNodesPerElement(type);
    //    type_p1 = Mesh::getP1ElementType(type);
    //    nb_nodes_per_element_p1[nb_good_types] = Mesh::getNbNodesPerElement(type_p1);

    conn_val[nb_good_types] = mesh.getConnectivity(type, _not_ghost).values;
    nb_element[nb_good_types] = mesh.getConnectivity(type, _not_ghost).getSize();
    nb_good_types++;
  }

  AKANTU_DEBUG_ASSERT(nb_good_types  != 0,
		      "Some elements must be found in right dimension to compute facets!");

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  //  node_offset.resize(nb_nodes + 1);
  //  UInt * node_offset_val = node_offset.values;

  /// count number of occurrence of each node
  //  memset(node_offset_val, 0, (nb_nodes + 1)*sizeof(UInt));
  for (UInt t = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element[t]; ++n) {
	++node_to_elem.rowOffset(conn_val[t][el_offset + n]);
	//	node_offset_val[conn_val[t][el_offset + n]]++;
      }
    }
  }

  // /// convert the occurrence array in a csr one
  // for (UInt i = 1; i < nb_nodes; ++i) node_offset_val[i] += node_offset_val[i-1];
  // for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  // node_offset_val[0] = 0;

  node_to_elem.countToCSR();
  node_to_elem.resizeCols();
  node_to_elem.beginInsertions();
  /// rearrange element to get the node-element list
  // node_to_elem.resize(node_offset_val[nb_nodes]);
  // UInt * node_to_elem_val = node_to_elem.values;

  for (UInt t = 0, linearized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linearized_el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element[t]; ++n)
	node_to_elem.insertInRow(conn_val[t][el_offset + n], linearized_el);
      // node_to_elem_val[node_offset_val[conn_val[t][el_offset + n]]++] = linearized_el;
    }

  node_to_elem.endInsertions();
  // for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  // node_offset_val[0] = 0;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2ElementsByElementType(const Mesh & mesh,
						ElementType type,
						CSR<UInt> & node_to_elem) {
						// Vector<UInt> & node_offset,
						// Vector<UInt> & node_to_elem) {

  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_elements = mesh.getConnectivity(type, _not_ghost).getSize();

  UInt * conn_val = mesh.getConnectivity(type, _not_ghost).values;

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  // node_offset.resize(nb_nodes + 1);
  // UInt * node_offset_val = node_offset.values;
  // memset(node_offset_val, 0, (nb_nodes + 1)*sizeof(UInt));

  /// count number of occurrence of each node
  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el*nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n)
      ++node_to_elem.rowOffset(conn_val[el_offset + n]);
      //      node_offset_val[conn_val[nb_nodes_per_element*el + n]]++;
  }

  /// convert the occurrence array in a csr one
  // for (UInt i = 1; i < nb_nodes; ++i) node_offset_val[i] += node_offset_val[i-1];
  // for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  // node_offset_val[0] = 0;
  node_to_elem.countToCSR();

  node_to_elem.resizeCols();
  node_to_elem.beginInsertions();

  /// save the element index in the node-element list
  // node_to_elem.resize(node_offset_val[nb_nodes]);
  // UInt * node_to_elem_val = node_to_elem.values;

  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el*nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      node_to_elem.insertInRow(conn_val[el_offset + n], el);
    }
    // node_to_elem_val[node_offset_val[conn_val[nb_nodes_per_element*el + n]]] = el;
    // node_offset_val[conn_val[nb_nodes_per_element*el + n]]++;
  }

  // ///  rearrange node_offset to start with 0
  // for (UInt i = nb_nodes; i > 0; --i) node_offset_val[i]  = node_offset_val[i-1];
  // node_offset_val[0] = 0;
  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacets(Mesh & mesh){
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  ByElementTypeReal barycenter;

  buildFacetsDimension(mesh, mesh, true, spatial_dimension, barycenter);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacets(Mesh & mesh, Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  ByElementTypeReal barycenter;

  /// generate facets
  buildFacetsDimension(mesh, mesh_facets, false, spatial_dimension, barycenter);

  /// compute their barycenters
  mesh_facets.initByElementTypeVector(barycenter, spatial_dimension, spatial_dimension - 1);
  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1);
  for(; it != end; ++it) {
    UInt nb_element = mesh_facets.getNbElement(*it);
    barycenter(*it).resize(nb_element);

    Vector<Real>::iterator<types::RVector> bary = barycenter(*it).begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> bary_end
      = barycenter(*it).end(spatial_dimension);

    for (UInt el = 0; bary != bary_end; ++bary, ++el) {
      mesh_facets.getBarycenter(el, *it, bary->storage());
    }
  }

  /// sort facets and generate subfacets
  for (UInt i = spatial_dimension - 1; i > 0; --i) {
    buildFacetsDimension(mesh_facets, mesh_facets, false, i, barycenter);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacetsDimension(Mesh & mesh,
				     Mesh & mesh_facets,
				     bool boundary_only,
				     UInt dimension,
				     ByElementTypeReal & barycenter){
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  UInt nb_good_types = 0;
  UInt spatial_dimension = mesh.getSpatialDimension();

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_facet[nb_types];

  UInt nb_facets[nb_types];
  UInt ** node_in_facet[nb_types];
  Vector<UInt> * connectivity_facets[nb_types];
  Vector<Vector<Element> > * element_to_subelement[nb_types];
  Vector<Element> * subelement_to_element[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];

  ElementType facet_type;

  Real epsilon = std::numeric_limits<Real>::epsilon();

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

    if(mesh.getSpatialDimension(type) != dimension) continue;

    nb_nodes_per_element[nb_good_types] = mesh.getNbNodesPerElement(type);

    facet_type = mesh.getFacetElementType(type);
    nb_facets[nb_good_types] = mesh.getNbFacetsPerElement(type);
    node_in_facet[nb_good_types] = mesh.getFacetLocalConnectivity(type);

    nb_nodes_per_facet[nb_good_types] = mesh.getNbNodesPerElement(facet_type);

    // getting connectivity of boundary facets
    connectivity_facets[nb_good_types] = mesh_facets.getConnectivityPointer(facet_type);
    connectivity_facets[nb_good_types]->resize(0);

    element_to_subelement[nb_good_types] = mesh_facets.getElementToSubelementPointer(facet_type);
    element_to_subelement[nb_good_types]->resize(0);

    conn_val[nb_good_types] = mesh.getConnectivity(type, _not_ghost).values;
    nb_element[nb_good_types] = mesh.getConnectivity(type, _not_ghost).getSize();

    subelement_to_element[nb_good_types] = mesh_facets.getSubelementToElementPointer(type);
    subelement_to_element[nb_good_types]->resize(nb_element[nb_good_types]);

    nb_good_types++;
  }


  CSR<UInt> node_to_elem;
  // Vector<UInt> node_offset;
  // Vector<UInt> node_to_elem;

  //  buildNode2Elements(mesh,node_offset,node_to_elem);
  buildNode2Elements(mesh, node_to_elem, dimension);

  // std::cout << "node offset " << std::endl << node_offset << std::endl;
  // std::cout << "node to elem " << std::endl << node_to_elem << std::endl;

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
	CSR<UInt>::iterator first_node_elements;
	//UInt * first_node_elements = node_to_elem.values+node_offset.values[facet_nodes[0]];
	UInt first_node_nelements = node_to_elem.getNbCols(facet_nodes[0]);
	  // node_offset.values[facet_nodes[0]+1]-
	  // node_offset.values[facet_nodes[0]];
	counter.resize(first_node_nelements);
	counter.clear();
	//loop over the other nodes to search intersecting elements,
	//which are the elements that share another node with the
	//starting element after first_node
	CSR<UInt>::iterator first_node_elements_end;
	first_node_elements_end = node_to_elem.end(facet_nodes[0]);
	for (UInt n = 1; n < nb_nodes_per_facet[t]; ++n) {
	  CSR<UInt>::iterator node_elements, node_elements_begin, node_elements_end;
	  node_elements_begin = node_to_elem.begin(facet_nodes[n]);
	  node_elements_end = node_to_elem.end(facet_nodes[n]);
	  //	  UInt * node_elements = node_to_elem.values+node_offset.values[facet_nodes[n]];
	    // node_offset.values[facet_nodes[n]+1]-
	    // node_offset.values[facet_nodes[n]];
	  UInt local_el = 0;
	  for (first_node_elements = node_to_elem.begin(facet_nodes[0]);
	       first_node_elements != first_node_elements_end;
	       ++first_node_elements, ++local_el) {
	    for (node_elements = node_elements_begin;
		 node_elements != node_elements_end;
		 ++node_elements) {
	      if (*first_node_elements == *node_elements) {
		++counter.values[local_el];
		// it may cause trouble:
		break;
	      }
	    }
	    //	    if (counter.values[local_el] == nb_nodes_per_facet[t]) break;
	  }
	}
	// bool connected_facet = false;

	//the connected elements are those for which counter is n_facet
	//	UInt connected_element = -1;

	// counting the number of elements connected to the facets and
	// taking the minimum element number, because the facet should
	// be inserted just once
	UInt nb_element_connected_to_facet = 0;
	UInt minimum_el_index = std::numeric_limits<UInt>::max();
	Vector<UInt> connected_elements;
	for (UInt el1 = 0; el1 < counter.getSize(); el1++) {
	  UInt el_index = node_to_elem(facet_nodes[0], el1);

	  if (counter.values[el1] == nb_nodes_per_facet[t]-1) {
	    ++nb_element_connected_to_facet;
	    minimum_el_index = std::min(minimum_el_index, el_index);
	    connected_elements.push_back(el_index);
	  }
	}

	if (minimum_el_index == linearized_el) {
	  if (!boundary_only || (boundary_only && nb_element_connected_to_facet == 1)) {
	    connectivity_facets[t]->push_back(facet_nodes);

	    UInt current_nb_facets = element_to_subelement[t]->getSize();
	    element_to_subelement[t]->resize(current_nb_facets + 1);
	    Vector<Element> & elements = (*element_to_subelement[t])(current_nb_facets);

	    // build elements_on_facets: linearized_el must come first
	    // in order to store the facet in the correct direction
	    // and avoid to invert the sign in the normal computation
	    elements.push_back(mesh.linearizedToElement(linearized_el));

	    /// boundary facet
	    if (nb_element_connected_to_facet == 1)
	      elements.push_back(ElementNull);
	    /// internal facet
	    else if (nb_element_connected_to_facet == 2)
	      elements.push_back(mesh.linearizedToElement(connected_elements.values[1]));
	    /// facet of facet
	    else {
	      UInt nb_connected = connected_elements.getSize();
	      for (UInt i = 1; i < nb_connected; ++i) {
		elements.push_back(mesh.linearizedToElement(connected_elements.values[i]));
	      }

	      /// check if sorting is needed:
	      /// - in 3D to sort triangles around segments
	      /// - in 2D to sort segments around points
	      if (dimension == spatial_dimension - 1) {

		/// barycentrical coordinates for each connected
		/// element with respect to start_node
		Vector<Real> connected_nodes(nb_connected, 2);

		const Vector<Real> & coord = mesh_facets.getNodes();
		const Vector<UInt> & facet_conn = mesh_facets.getConnectivity(facet_type);

		/// node around which the sorting is carried out is
		/// the first node of the current facet
		UInt start_node = facet_conn(facet_conn.getSize()-1, 0);
		Real start_coord[spatial_dimension];
		for (UInt dim = 0; dim < spatial_dimension; ++dim) {
		  start_coord[dim] = coord(start_node, dim);
		}

		if (spatial_dimension == 3) {
		  /// vector connecting facet first node to second
		  Real tangent[spatial_dimension];
		  /// vector connecting facet first node and
		  /// barycenter of elements(0)
		  Real temp[spatial_dimension];
		  /// two normals to the segment facet to define the
		  /// reference system
		  Real normal1[spatial_dimension];
		  Real normal2[spatial_dimension];

		  Vector<Real> & bar = barycenter(elements(0).type);

		  /// facet second node
		  UInt second_node = facet_conn(facet_conn.getSize()-1, 1);

		  /// construction of tangent and temp arrays
		  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
		    Real x1, x2;
		    x1 = coord(second_node, dim);
		    x2 = bar(elements(0).element, dim);
		    tangent[dim] = x1 - start_coord[dim];
		    temp[dim] = x2 - start_coord[dim];
		  }

		  /// get normal1 and normal2
		  Math::normalize3(tangent);
		  Math::vectorProduct3(tangent, temp, normal1);
		  Math::normalize3(normal1);
		  Math::vectorProduct3(tangent, normal1, normal2);

		  for (UInt n = 0; n < nb_connected; ++n) {
		    Real bary_coord[spatial_dimension];
		    Vector<Real> & bary = barycenter(elements(n).type);

		    /// get the barycenter local coordinates
		    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
		      bary_coord[dim] = bary(elements(n).element, dim)
			- start_coord[dim];
		    }
		    /// project the barycenter coordinates on the two
		    /// normals to have them on the same plane
		    connected_nodes(n, 0) = Math::vectorDot(bary_coord, normal1, spatial_dimension);
		    connected_nodes(n, 1) = Math::vectorDot(bary_coord, normal2, spatial_dimension);
		  }
		}
		else if (spatial_dimension == 2) {
		  for (UInt n = 0; n < nb_connected; ++n) {
		    Vector<Real> & bary = barycenter(elements(n).type);
		    /// get the barycenter local coordinates
		    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
		      connected_nodes(n, dim) = bary(elements(n).element, dim)
			- start_coord[dim];
		    }
		  }
		}

		/// associate to each element a real value based on
		/// atan2 function (check wikipedia)
		std::map<Element, Real, CompElementLess> atan2;

		for (UInt n = 0; n < nb_connected; ++n) {
		  Real x = connected_nodes(n, 0);
		  Real y = connected_nodes(n, 1);
		  /// in order to avoid division by zero:
		  if (std::abs(y) <= std::abs(y) * epsilon && x < 0)
		    y = Math::getTolerance();
		  atan2[elements(n)] = y / (sqrt(x * x + y * y) + x);
		}

		/// sort elements according to their atan2 values
		ElementSorter sorter(atan2);
		std::sort(elements.storage(), elements.storage() + elements.getSize(), sorter);

	      }
	    }

	    /// current facet index
	    UInt current_facet = connectivity_facets[t]->getSize() - 1;

	    /// loop on every element connected to current facet and
	    /// insert current facet in the first free spot of the
	    /// subelement_to_element vector
	    for (UInt elem = 0; elem < elements.getSize(); ++elem) {
	      if (elements(elem).type != _not_defined) {
		for (UInt f_in = 0; f_in < nb_facets[t]; ++f_in) {
		  if ((*subelement_to_element[t])(elements(elem).element, f_in).type == _not_defined) {
		    (*subelement_to_element[t])(elements(elem).element, f_in).type = facet_type;
		    (*subelement_to_element[t])(elements(elem).element, f_in).element = current_facet;
		    break;
		  }
		}
	      }
	    }

	  }
	}

      }
    }
  }

  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
// void MeshUtils::buildNormals(Mesh & mesh,UInt spatial_dimension){
//   AKANTU_DEBUG_IN();
//   const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
//   Mesh::ConnectivityTypeList::const_iterator it;

//   UInt nb_types = type_list.size();
//   UInt nb_nodes_per_element[nb_types];
//   UInt nb_nodes_per_element_p1[nb_types];

//   UInt nb_good_types = 0;

//   Vector<UInt> * connectivity[nb_types];
//   Vector<Real> * normals[nb_types];

//   if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();

//   for(it = type_list.begin(); it != type_list.end(); ++it) {
//     ElementType type = *it;
//     ElementType type_p1 = mesh.getP1ElementType(type);
//     if(mesh.getSpatialDimension(type) != spatial_dimension) continue;

//     nb_nodes_per_element[nb_good_types]    = mesh.getNbNodesPerElement(type);
//     nb_nodes_per_element_p1[nb_good_types] = mesh.getNbNodesPerElement(type_p1);

//     // getting connectivity
//     connectivity[nb_good_types] = mesh.getConnectivityPointer(type);
//     if (!connectivity[nb_good_types])
//       AKANTU_DEBUG_ERROR("connectivity is not allocatted : this should probably not have happened");

//     //getting array of normals
//     normals[nb_good_types] = mesh.getNormalsPointer(type);
//     if(normals[nb_good_types])
//       normals[nb_good_types]->resize(0);
//     else
//       normals[nb_good_types] = &mesh.createNormals(type);

//     nb_good_types++;
//   }
//   AKANTU_DEBUG_OUT();
// }

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

  std::map<UInt, UInt> renumbering_map;

  /// renumber the nodes
  for (UInt el = 0; el < nb_local_element + nb_ghost_element; ++el) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      std::map<UInt, UInt>::iterator it = renumbering_map.find(*conn);
      //      Int nn = nodes_numbers.find(*conn);
      if(it == renumbering_map.end()) {
	UInt old_node = *conn;
	nodes_numbers.push_back(old_node);
	*conn = nodes_numbers.getSize() - 1;

	renumbering_map[old_node] = *conn;
      } else {
	*conn = it->second;
      }
      conn++;
    }
  }

  renumbering_map.clear();

  /// copy the renumbered connectivity to the right place
  Vector<UInt> * local_conn = mesh.getConnectivityPointer(type);
  local_conn->resize(nb_local_element);
  memcpy(local_conn->values,
	 local_connectivities,
	 nb_local_element * nb_nodes_per_element * sizeof(UInt));

  Vector<UInt> * ghost_conn = mesh.getConnectivityPointer(type,_ghost);
  ghost_conn->resize(nb_ghost_element);
  memcpy(ghost_conn->values,
	 local_connectivities + nb_local_element * nb_nodes_per_element,
	 nb_ghost_element * nb_nodes_per_element * sizeof(UInt));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::setUIntData(Mesh & mesh, UInt * data, UInt nb_tags, const ElementType & type) {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type, _not_ghost);
  UInt nb_ghost_element = mesh.getNbElement(type, _ghost);

  char * names = reinterpret_cast<char *>(data + (nb_element + nb_ghost_element) * nb_tags);
  UIntDataMap & uint_data_map = mesh.getUIntDataMap(type, _not_ghost);
  UIntDataMap & ghost_uint_data_map = mesh.getUIntDataMap(type, _ghost);

  for (UInt t = 0; t < nb_tags; ++t) {
    std::string name(names);
    //    std::cout << name << std::endl;
    names += name.size() + 1;

    UIntDataMap::iterator it = uint_data_map.find(name);
    if(it == uint_data_map.end()) {
      uint_data_map[name] = new Vector<UInt>(0, 1, name);
      it = uint_data_map.find(name);
    }
    it->second->resize(nb_element);
    memcpy(it->second->values, data, nb_element * sizeof(UInt));
    data += nb_element;

    it = ghost_uint_data_map.find(name);
    if(it == ghost_uint_data_map.end()) {
      ghost_uint_data_map[name] = new Vector<UInt>(0, 1, name);
      it = ghost_uint_data_map.find(name);
    }
    it->second->resize(nb_ghost_element);
    memcpy(it->second->values, data, nb_ghost_element * sizeof(UInt));
    data += nb_ghost_element;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildSurfaceID(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  // Vector<UInt> node_offset;
  // Vector<UInt> node_to_elem;

  CSR<UInt> node_to_elem;

  /// Get list of surface elements
  UInt spatial_dimension = mesh.getSpatialDimension();

  //  buildNode2Elements(mesh, node_offset, node_to_elem, spatial_dimension-1);
  buildNode2Elements(mesh, node_to_elem, spatial_dimension-1);

  /// Find which types of elements have been linearized
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  ElementType lin_element_type[nb_types];
  UInt nb_lin_types = 0;

  UInt nb_nodes_per_element[nb_types];
  UInt nb_nodes_per_element_p1[nb_types];

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types+1];

  ElementType type_p1;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    ElementType facet_type = mesh.getFacetElementType(type);
    lin_element_type[nb_lin_types] = facet_type;
    nb_nodes_per_element[nb_lin_types]    = Mesh::getNbNodesPerElement(facet_type);
    type_p1 = Mesh::getP1ElementType(facet_type);
    nb_nodes_per_element_p1[nb_lin_types] = Mesh::getNbNodesPerElement(type_p1);

    conn_val[nb_lin_types] = mesh.getConnectivity(facet_type, _not_ghost).values;
    nb_element[nb_lin_types] = mesh.getNbElement(facet_type, _not_ghost);
    nb_lin_types++;
  }

  for (UInt i = 1; i < nb_lin_types; ++i) nb_element[i] += nb_element[i+1];
  for (UInt i = nb_lin_types; i > 0; --i) nb_element[i] = nb_element[i-1];
  nb_element[0] = 0;

  /// Find close surfaces
  Vector<Int> surface_value_id(1, nb_element[nb_lin_types], -1);
  Int * surf_val = surface_value_id.values;
  UInt nb_surfaces = 0;

  UInt nb_cecked_elements;
  UInt nb_elements_to_ceck;
  UInt * elements_to_ceck = new UInt [nb_element[nb_lin_types]];
  memset(elements_to_ceck, 0, nb_element[nb_lin_types]*sizeof(UInt));

  for (UInt lin_el = 0; lin_el < nb_element[nb_lin_types]; ++lin_el) {

    if(surf_val[lin_el] != -1) continue; /* Surface id already assigned */

    /* First element of new surface */
    surf_val[lin_el] = nb_surfaces;
    nb_cecked_elements = 0;
    nb_elements_to_ceck = 1;
    memset(elements_to_ceck, 0, nb_element[nb_lin_types]*sizeof(UInt));
    elements_to_ceck[0] = lin_el;

    // Find others elements belonging to this surface
    while(nb_cecked_elements < nb_elements_to_ceck) {

      UInt ceck_lin_el = elements_to_ceck[nb_cecked_elements];

      // Transform linearized index of element into ElementType one
      UInt lin_type_nb = 0;
      while (ceck_lin_el >= nb_element[lin_type_nb+1])
	lin_type_nb++;
      UInt ceck_el = ceck_lin_el - nb_element[lin_type_nb];

      // Get connected elements
      UInt el_offset = ceck_el*nb_nodes_per_element[lin_type_nb];
      for (UInt n = 0; n < nb_nodes_per_element_p1[lin_type_nb]; ++n) {
	UInt node_id = conn_val[lin_type_nb][el_offset + n];
	CSR<UInt>::iterator it_n;
	for (it_n = node_to_elem.begin(node_id); it_n != node_to_elem.end(node_id); ++it_n) {
	  //	for (UInt i = node_offset.values[node_id]; i < node_offset.values[node_id+1]; ++i) {
	  if(surf_val[*it_n] == -1) { /* Found new surface element */
	    surf_val[*it_n] = nb_surfaces;
	    elements_to_ceck[nb_elements_to_ceck] = *it_n;
	    nb_elements_to_ceck++;
	  }
	  // if(surf_val[node_to_elem.values[i]] == -1) { /* Found new surface element */
	  //   surf_val[node_to_elem.values[i]] = nb_surfaces;
	  //   elements_to_ceck[nb_elements_to_ceck] = node_to_elem.values[i];
	  //   nb_elements_to_ceck++;
	  // }
	}
      }

      nb_cecked_elements++;
    }

    nb_surfaces++;
  }

  delete [] elements_to_ceck;

  /// Transform local linearized element index in the global one
  for (UInt i = 0; i < nb_lin_types; ++i) nb_element[i] = nb_element[i+1] - nb_element[i];
  UInt el_offset = 0;

  for (UInt type_it = 0; type_it < nb_lin_types; ++type_it) {
    ElementType type = lin_element_type[type_it];
    Vector<UInt> * surf_id_type = mesh.getSurfaceIDPointer(type, _not_ghost);
    surf_id_type->resize(nb_element[type_it]);
    surf_id_type->clear();
    for (UInt el = 0; el < nb_element[type_it]; ++el)
      surf_id_type->values[el] = surf_val[el+el_offset];
    el_offset += nb_element[type_it];
  }

  /// Set nb_surfaces in mesh
  mesh.nb_surfaces = nb_surfaces;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNodesPerSurface(const Mesh & mesh, CSR<UInt> & nodes_per_surface) {
  AKANTU_DEBUG_IN();

  UInt nb_surfaces = mesh.getNbSurfaces();
  UInt nb_nodes    = mesh.getNbNodes();
  UInt spatial_dimension = mesh.getSpatialDimension(); //surface elements

  UInt  nb_facet_types = 0;
  ElementType facet_type[_max_element_type];

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);
  for(; it != end; ++it) {
    facet_type[nb_facet_types++] = mesh.getFacetElementType(*it);
  }

  UInt * surface_nodes_id = new UInt[nb_nodes*nb_surfaces];
  std::fill_n(surface_nodes_id, nb_surfaces*nb_nodes, 0);

  for(UInt t = 0; t < nb_facet_types; ++t) {
    ElementType type = facet_type[t];

    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    UInt * connecticity = mesh.getConnectivity(type, _not_ghost).values;
    UInt * surface_id = mesh.getSurfaceID(type, _not_ghost).values;;

    for (UInt el = 0; el < nb_element; ++el) {
      UInt offset = *surface_id * nb_nodes;
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	surface_nodes_id[offset + *connecticity] = 1;
	++connecticity;
      }
      ++surface_id;
    }
  }

  nodes_per_surface.resizeRows(nb_surfaces);
  nodes_per_surface.clearRows();

  UInt * surface_nodes_id_tmp = surface_nodes_id;
  for (UInt s = 0; s < nb_surfaces; ++s)
    for (UInt n = 0; n < nb_nodes; ++n)
      nodes_per_surface.rowOffset(s) += *surface_nodes_id_tmp++;

  nodes_per_surface.countToCSR();

  nodes_per_surface.resizeCols();
  nodes_per_surface.beginInsertions();
  surface_nodes_id_tmp = surface_nodes_id;
  for (UInt s = 0; s < nb_surfaces; ++s)
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (*surface_nodes_id_tmp == 1) nodes_per_surface.insertInRow(s, n);
      surface_nodes_id_tmp++;
    }
  nodes_per_surface.endInsertions();

  delete [] surface_nodes_id;

  AKANTU_DEBUG_OUT();
}



__END_AKANTU__

//  LocalWords:  ElementType
