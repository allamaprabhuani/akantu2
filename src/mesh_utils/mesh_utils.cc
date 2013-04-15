/**
 * @file   mesh_utils.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Fri Aug 20 12:19:44 2010
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
#include "aka_safe_enum.hh"
#include "fem.hh"


__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2Elements(const Mesh & mesh,
				   CSR<Element> & node_to_elem,
				   UInt spatial_dimension) {
  AKANTU_DEBUG_IN();
  if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();


  /// count number of occurrence of each node
  UInt nb_nodes = mesh.getNbNodes();

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  AKANTU_DEBUG_ASSERT(mesh.firstType(spatial_dimension) !=
		      mesh.lastType(spatial_dimension),
		      "Some elements must be found in right dimension to compute facets!");

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator first = mesh.firstType(spatial_dimension, *gt);
    Mesh::type_iterator last  = mesh.lastType(spatial_dimension, *gt);

    for (; first != last; ++first) {
      ElementType type = *first;
      UInt nb_element = mesh.getNbElement(type, *gt);
      Array<UInt>::const_iterator< Vector<UInt> > conn_it =
	mesh.getConnectivity(type, *gt).begin(Mesh::getNbNodesPerElement(type));

      for (UInt el = 0; el < nb_element; ++el, ++conn_it)
	for (UInt n = 0; n < conn_it->size(); ++n)
	  ++node_to_elem.rowOffset((*conn_it)(n));
    }
  }

  node_to_elem.countToCSR();
  node_to_elem.resizeCols();


  /// rearrange element to get the node-element list
  Element e;
  node_to_elem.beginInsertions();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    Mesh::type_iterator first = mesh.firstType(spatial_dimension, *gt);
    Mesh::type_iterator last  = mesh.lastType(spatial_dimension, *gt);
    e.ghost_type = *gt;
    for (; first != last; ++first) {
      ElementType type = *first;
      e.type = type;
      UInt nb_element = mesh.getNbElement(type, *gt);
      Array<UInt>::const_iterator< Vector<UInt> > conn_it =
	mesh.getConnectivity(type, *gt).begin(Mesh::getNbNodesPerElement(type));

      for (UInt el = 0; el < nb_element; ++el, ++conn_it) {
	e.element = el;
	for (UInt n = 0; n < conn_it->size(); ++n)
	  node_to_elem.insertInRow((*conn_it)(n), e);
      }
    }
  }

  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * This function should disappear in the future
 */
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

  UInt * conn_val[nb_types];
  UInt nb_element[nb_types];


  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(Mesh::getSpatialDimension(type) != spatial_dimension) continue;

    nb_nodes_per_element[nb_good_types]    = Mesh::getNbNodesPerElement(type);
    conn_val[nb_good_types] = mesh.getConnectivity(type, _not_ghost).values;
    nb_element[nb_good_types] = mesh.getConnectivity(type, _not_ghost).getSize();
    nb_good_types++;
  }

  AKANTU_DEBUG_ASSERT(nb_good_types  != 0,
		      "Some elements must be found in right dimension to compute facets!");

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  /// count number of occurrence of each node
  for (UInt t = 0; t < nb_good_types; ++t) {
    for (UInt el = 0; el < nb_element[t]; ++el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element[t]; ++n) {
	++node_to_elem.rowOffset(conn_val[t][el_offset + n]);
      }
    }
  }

  node_to_elem.countToCSR();
  node_to_elem.resizeCols();
  node_to_elem.beginInsertions();

  /// rearrange element to get the node-element list
  for (UInt t = 0, linearized_el = 0; t < nb_good_types; ++t)
    for (UInt el = 0; el < nb_element[t]; ++el, ++linearized_el) {
      UInt el_offset = el*nb_nodes_per_element[t];
      for (UInt n = 0; n < nb_nodes_per_element[t]; ++n)
	node_to_elem.insertInRow(conn_val[t][el_offset + n], linearized_el);
    }

  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildNode2ElementsByElementType(const Mesh & mesh,
						ElementType type,
						CSR<UInt> & node_to_elem) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_elements = mesh.getConnectivity(type, _not_ghost).getSize();

  UInt * conn_val = mesh.getConnectivity(type, _not_ghost).values;

  /// array for the node-element list
  node_to_elem.resizeRows(nb_nodes);
  node_to_elem.clearRows();

  /// count number of occurrence of each node
  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el*nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n)
      ++node_to_elem.rowOffset(conn_val[el_offset + n]);
  }

  /// convert the occurrence array in a csr one
  node_to_elem.countToCSR();

  node_to_elem.resizeCols();
  node_to_elem.beginInsertions();

  /// save the element index in the node-element list
  for (UInt el = 0; el < nb_elements; ++el) {
    UInt el_offset = el*nb_nodes_per_element;
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      node_to_elem.insertInRow(conn_val[el_offset + n], el);
    }
  }

  node_to_elem.endInsertions();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacets(Mesh & mesh){
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  ByElementTypeReal barycenter;
  ByElementTypeUInt prank_to_element;

  buildFacetsDimension(mesh,
		       mesh,
		       true,
		       spatial_dimension,
		       barycenter,
		       prank_to_element);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacets(Mesh & mesh,
			       Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  ByElementTypeUInt prank_to_element;
  buildAllFacetsParallel(mesh, mesh_facets, prank_to_element);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacetsParallel(Mesh & mesh,
				       Mesh & mesh_facets,
				       ByElementTypeUInt & prank_to_element) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  ByElementTypeReal barycenter;

  /// generate facets
  buildFacetsDimension(mesh,
		       mesh_facets,
		       false,
		       spatial_dimension,
		       barycenter,
		       prank_to_element);

  /// compute their barycenters
  mesh_facets.initByElementTypeArray(barycenter,
				     spatial_dimension,
				     spatial_dimension - 1);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, *gt);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, *gt);

    for(; it != end; ++it) {
      UInt nb_element = mesh_facets.getNbElement(*it, *gt);
      barycenter(*it, *gt).resize(nb_element);

      Array<Real>::iterator< Vector<Real> > bary
	= barycenter(*it, *gt).begin(spatial_dimension);
      Array<Real>::iterator< Vector<Real> > bary_end
	= barycenter(*it, *gt).end(spatial_dimension);

      for (UInt el = 0; bary != bary_end; ++bary, ++el) {
	mesh_facets.getBarycenter(el, *it, bary->storage(), *gt);
      }
    }
  }

  /// sort facets and generate subfacets
  for (UInt i = spatial_dimension - 1; i > 0; --i) {
    buildFacetsDimension(mesh_facets,
			 mesh_facets,
			 false,
			 i,
			 barycenter,
			 prank_to_element);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void MeshUtils::buildFacetsDimension(Mesh & mesh,
				     Mesh & mesh_facets,
				     bool boundary_only,
				     UInt dimension,
				     ByElementTypeReal & barycenter,
				     ByElementTypeUInt & prank_to_element){
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  const Array<Real> & mesh_facets_nodes = mesh_facets.getNodes();
  const Array<Real>::const_iterator< Vector<Real> > mesh_facets_nodes_it =
    mesh_facets_nodes.begin(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    Mesh::type_iterator first = mesh.firstType(dimension, ghost_type);
    Mesh::type_iterator last  = mesh.lastType(dimension, ghost_type);
    for(; first != last; ++first) {
      ElementType type = *first;
      if(mesh.getSpatialDimension(type) != dimension) continue;

      ElementType facet_type = mesh.getFacetType(type);
      UInt nb_element = mesh.getNbElement(type, ghost_type);

      // getting connectivity of boundary facets
      Array<UInt> * connectivity_facets = mesh_facets.getConnectivityPointer(facet_type, ghost_type);
      connectivity_facets->resize(0);
      Array< std::vector<Element> > * element_to_subelement =
	mesh_facets.getElementToSubelementPointer(facet_type, ghost_type);
      element_to_subelement->resize(0);
      Array<Element> * subelement_to_element = mesh_facets.getSubelementToElementPointer(type, ghost_type);
      subelement_to_element->resize(nb_element);
    }
  }

  CSR<Element> node_to_elem;
  buildNode2Elements(mesh, node_to_elem, dimension);

  Array<UInt> counter;

  const Array<Int> * nodes_type = NULL;
  nodes_type = &(mesh.getNodesType());

  Element current_element;
  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    GhostType facet_ghost_type = ghost_type;

    current_element.ghost_type = ghost_type;
    Mesh::type_iterator first = mesh.firstType(dimension, ghost_type);
    Mesh::type_iterator last  = mesh.lastType(dimension, ghost_type);

    for(; first != last; ++first) {
      ElementType type = *first;
      ElementType facet_type = mesh.getFacetType(type);

      current_element.type = type;

      UInt nb_element = mesh.getNbElement(type, ghost_type);
      Array< std::vector<Element> > * element_to_subelement = &mesh_facets.getElementToSubelement(facet_type, ghost_type);
      Array<UInt> * connectivity_facets = &mesh_facets.getConnectivity(facet_type, ghost_type);

      for (UInt el = 0; el < nb_element; ++el) {
	current_element.element = el;
	Matrix<UInt> facets = mesh.getFacetConnectivity(el, type, ghost_type);
	UInt nb_nodes_per_facet = facets.cols();

	for (UInt f = 0; f < facets.rows(); ++f) {
	  Vector<UInt> facet(nb_nodes_per_facet);
	  for (UInt n = 0; n < nb_nodes_per_facet; ++n) facet(n) = facets(f, n);

	  UInt first_node_nb_elements = node_to_elem.getNbCols(facets(f, 0));
	  counter.resize(first_node_nb_elements);
	  counter.clear();

	  //loop over the other nodes to search intersecting elements,
	  //which are the elements that share another node with the
	  //starting element after first_node
	  CSR<Element>::iterator first_node_elements = node_to_elem.begin(facet(0));
	  CSR<Element>::iterator first_node_elements_end = node_to_elem.end(facet(0));
	  UInt local_el = 0;
	  for (; first_node_elements != first_node_elements_end; ++first_node_elements, ++local_el) {
	    for (UInt n = 1; n < nb_nodes_per_facet; ++n) {
	      CSR<Element>::iterator node_elements_begin = node_to_elem.begin(facet(n));
	      CSR<Element>::iterator node_elements_end   = node_to_elem.end  (facet(n));
	      counter(local_el) += std::count(node_elements_begin, node_elements_end, *first_node_elements);
	    }
	  }

	  // counting the number of elements connected to the facets and
	  // taking the minimum element number, because the facet should
	  // be inserted just once
	  UInt nb_element_connected_to_facet = 0;
	  Element minimum_el = ElementNull;
	  Array<Element> connected_elements;
	  for (UInt el_f = 0; el_f < first_node_nb_elements; el_f++) {
	    Element real_el = node_to_elem(facet(0), el_f);
	    if (counter(el_f) == nb_nodes_per_facet - 1) {
	      ++nb_element_connected_to_facet;
	      minimum_el = std::min(minimum_el, real_el);
	      connected_elements.push_back(real_el);
	    }
	  }

	  if (minimum_el == current_element) {
	    bool full_ghost_facet = false;

	    if (nodes_type != NULL) {
	      UInt n = 0;
	      while (n < nb_nodes_per_facet && (*nodes_type)(facet(n)) == -3) ++n;
	      if (n == nb_nodes_per_facet) full_ghost_facet = true;
	    }

	    if (!full_ghost_facet) {
	      if (!boundary_only || (boundary_only && nb_element_connected_to_facet == 1)) {

		std::vector<Element> elements;

		// build elements_on_facets: linearized_el must come first
		// in order to store the facet in the correct direction
		// and avoid to invert the sign in the normal computation
		elements.push_back(current_element);

		/// boundary facet
		if (nb_element_connected_to_facet == 1)
		  elements.push_back(ElementNull);
		/// internal facet
		else if (nb_element_connected_to_facet == 2) {
		  elements.push_back(connected_elements(1));

		  /// check if facet is in between ghost and normal
		  /// elements: if it's the case, the facet is either
		  /// ghost or not ghost. The criterion to decide this
		  /// is arbitrary. It was chosen to check the processor
		  /// id (prank) of the two neighboring elements. If
		  /// prank of the ghost element is lower than prank of
		  /// the normal one, the facet is not ghost, otherwise
		  /// it's ghost
		  GhostType gt[2];
		  for (UInt el = 0; el < connected_elements.getSize(); ++el)
		    gt[el] = connected_elements(el).ghost_type;

		  if (gt[0] + gt[1] == 1) {
		    try {
		      UInt prank[2];
		      for (UInt el = 0; el < 2; ++el) {
			UInt current_el = connected_elements(el).element;
			ElementType current_type = connected_elements(el).type;
			GhostType current_gt = connected_elements(el).ghost_type;

			Array<UInt> & prank_to_el = prank_to_element(current_type,
								     current_gt);

			prank[el] = prank_to_el(current_el);
		      }

		      bool ghost_one = (gt[0] != _ghost);

		      if (prank[ghost_one] > prank[!ghost_one])
			facet_ghost_type = _not_ghost;
		      else
			facet_ghost_type = _ghost;

		      connectivity_facets = &mesh_facets.getConnectivity(facet_type, facet_ghost_type);
		      element_to_subelement = &mesh_facets.getElementToSubelement(facet_type, facet_ghost_type);
		    } catch (...) { }
		  }
		}
		/// facet of facet
		else {
		  for (UInt i = 1; i < nb_element_connected_to_facet; ++i) {
		    elements.push_back(connected_elements(i));
		  }

		  /// check if sorting is needed:
		  /// - in 3D to sort triangles around segments
		  /// - in 2D to sort segments around points
		  if (dimension == spatial_dimension - 1) {
		    /// node around which the sorting is carried out is
		    /// the first node of the current facet
		    const Vector<Real> & first_node_coord = mesh_facets_nodes_it[facet(0)];

		    /// associate to each element a real value based on
		    /// atan2 function (check wikipedia)
		    std::map<Element, Real, CompElementLess> atan2;

		    if (spatial_dimension == 3) {
		      const Vector<Real> & second_node_coord = mesh_facets_nodes_it[facet(1)];

		      /// vector connecting facet first node to second
		      Vector<Real> tangent(spatial_dimension);
		      tangent = second_node_coord;
		      tangent -= first_node_coord;
		      tangent.normalize();

		      const Array<Real>::const_iterator< Vector<Real> > bar =
			barycenter(elements[0].type, elements[0].ghost_type).begin(spatial_dimension);

		      /// vector connecting facet first node and
		      /// barycenter of elements(0)
		      Vector<Real> bary_coord(spatial_dimension);
		      bary_coord.copy(bar[elements[0].element]);
		      bary_coord -= first_node_coord;

		      /// two normals to the segment facet to define the
		      /// reference system
		      Vector<Real> normal1(spatial_dimension);
		      Vector<Real> normal2(spatial_dimension);

		      /// get normal1 and normal2
		      normal1.crossProduct(tangent, bary_coord);
		      normal1.normalize();
		      normal2.crossProduct(tangent, normal1);

		      /// project the barycenter coordinates on the two
		      /// normals to have them on the same plane
		      atan2[elements[0]] = std::atan2(bary_coord.dot(normal2), bary_coord.dot(normal1));

		      for (UInt n = 1; n < nb_element_connected_to_facet; ++n) {
			const Array<Real>::const_iterator< Vector<Real> > bar_it =
			  barycenter(elements[n].type, elements[n].ghost_type).begin(spatial_dimension);
			bary_coord.copy(bar_it[elements[n].element]);
			bary_coord -= first_node_coord;

			/// project the barycenter coordinates on the two
			/// normals to have them on the same plane
			atan2[elements[n]] = std::atan2(bary_coord.dot(normal2), bary_coord.dot(normal1));
		      }
		    }
		    else if (spatial_dimension == 2) {
		      for (UInt n = 0; n < nb_element_connected_to_facet; ++n) {
			const Array<Real>::const_iterator< Vector<Real> > bar_it =
			  barycenter(elements[n].type, elements[n].ghost_type).begin(spatial_dimension);
			Vector<Real> bary_coord(spatial_dimension);
			bary_coord.copy(bar_it[elements[n].element]);
			bary_coord -= first_node_coord;
			atan2[elements[n]] = std::atan2(bary_coord(1), bary_coord(0));
		      }
		    }

		    /// sort elements according to their atan2 values
		    ElementSorter sorter(atan2);
		    std::sort(elements.begin(), elements.end(), sorter);
		  }
		}

		element_to_subelement->push_back(elements);
		connectivity_facets->push_back(facet);

		/// current facet index
		UInt current_facet = connectivity_facets->getSize() - 1;

		/// loop on every element connected to current facet and
		/// insert current facet in the first free spot of the
		/// subelement_to_element vector
		for (UInt elem = 0; elem < elements.size(); ++elem) {
		  Element loc_el = elements[elem];

		  if (loc_el.type != _not_defined) {
		    Array<Element> & subelement_to_element =
		      mesh_facets.getSubelementToElement(type, loc_el.ghost_type);

		    for (UInt f_in = 0; f_in < facets.rows(); ++f_in) {
		      if (subelement_to_element(loc_el.element, f_in).type == _not_defined) {
			subelement_to_element(loc_el.element, f_in).type = facet_type;
			subelement_to_element(loc_el.element, f_in).element = current_facet;
			subelement_to_element(loc_el.element, f_in).ghost_type = facet_ghost_type;
			break;
		      }
		    }
		  }
		}

		/// reset connectivity in case a facet was found in
		/// between ghost and normal elements
		if (facet_ghost_type != ghost_type) {
		  facet_ghost_type = ghost_type;
		  connectivity_facets = &mesh_facets.getConnectivity(facet_type, facet_ghost_type);
		  element_to_subelement = &mesh_facets.getElementToSubelement(facet_type, facet_ghost_type);
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
void MeshUtils::renumberMeshNodes(Mesh & mesh,
				  UInt * local_connectivities,
				  UInt nb_local_element,
				  UInt nb_ghost_element,
				  ElementType type,
				  Array<UInt> & old_nodes_numbers) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  std::map<UInt, UInt> renumbering_map;
  for (UInt i = 0; i < old_nodes_numbers.getSize(); ++i) {
    renumbering_map[old_nodes_numbers(i)] = i;
  }

  /// renumber the nodes
  renumberNodesInConnectivity(local_connectivities,
			      (nb_local_element + nb_ghost_element)*nb_nodes_per_element,
			      renumbering_map);

  std::map<UInt, UInt>::iterator it = renumbering_map.begin();
  std::map<UInt, UInt>::iterator end = renumbering_map.end();
  old_nodes_numbers.resize(renumbering_map.size());
  for (;it != end; ++it) {
    old_nodes_numbers(it->second) = it->first;
  }
  renumbering_map.clear();

  /// copy the renumbered connectivity to the right place
  Array<UInt> * local_conn = mesh.getConnectivityPointer(type);
  local_conn->resize(nb_local_element);
  memcpy(local_conn->values,
	 local_connectivities,
	 nb_local_element * nb_nodes_per_element * sizeof(UInt));

  Array<UInt> * ghost_conn = mesh.getConnectivityPointer(type, _ghost);
  ghost_conn->resize(nb_ghost_element);
  memcpy(ghost_conn->values,
	 local_connectivities + nb_local_element * nb_nodes_per_element,
	 nb_ghost_element * nb_nodes_per_element * sizeof(UInt));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::renumberNodesInConnectivity(UInt * list_nodes,
					    UInt nb_nodes,
					    std::map<UInt, UInt> & renumbering_map) {
  AKANTU_DEBUG_IN();

  UInt * connectivity = list_nodes;
  UInt new_node_num = renumbering_map.size();
  for (UInt n = 0; n < nb_nodes; ++n, ++connectivity) {
    UInt & node = *connectivity;
    std::map<UInt, UInt>::iterator it = renumbering_map.find(node);
    if(it == renumbering_map.end()) {
      UInt old_node = node;
      renumbering_map[old_node] = new_node_num;
      node = new_node_num;
      ++new_node_num;
    } else {
      node = it->second;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::purifyMesh(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  std::map<UInt, UInt> renumbering_map;

  RemovedNodesEvent remove_nodes(mesh);
  Array<UInt> & nodes_removed = remove_nodes.getList();

  for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;

    Mesh::type_iterator it  = mesh.firstType(0, ghost_type, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(0, ghost_type, _ek_not_defined);
    for(; it != end; ++it) {

      ElementType type(*it);
      UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

      const Array<UInt> & connectivity_vect = mesh.getConnectivity(type, ghost_type);
      UInt nb_element(connectivity_vect.getSize());
      UInt * connectivity = connectivity_vect.storage();

      renumberNodesInConnectivity (connectivity, nb_element*nb_nodes_per_element, renumbering_map);
    }
  }

  Array<UInt> & new_numbering = remove_nodes.getNewNumbering();
  std::fill(new_numbering.begin(), new_numbering.end(), UInt(-1));

  std::map<UInt, UInt>::iterator it = renumbering_map.begin();
  std::map<UInt, UInt>::iterator end = renumbering_map.end();
  for (; it != end; ++it) {
    new_numbering(it->first) = it->second;
  }


  for (UInt i = 0; i < new_numbering.getSize(); ++i) {
    if(new_numbering(i) == UInt(-1))
      nodes_removed.push_back(i);
  }

  mesh.sendEvent(remove_nodes);

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
      uint_data_map[name] = new Array<UInt>(0, 1, name);
      it = uint_data_map.find(name);
    }
    it->second->resize(nb_element);
    memcpy(it->second->values, data, nb_element * sizeof(UInt));
    data += nb_element;

    it = ghost_uint_data_map.find(name);
    if(it == ghost_uint_data_map.end()) {
      ghost_uint_data_map[name] = new Array<UInt>(0, 1, name);
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

  // Array<UInt> node_offset;
  // Array<UInt> node_to_elem;

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

    ElementType facet_type = mesh.getFacetType(type);
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
  Array<Int> surface_value_id(1, nb_element[nb_lin_types], -1);
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
    Array<UInt> * surf_id_type = mesh.getSurfaceIDPointer(type, _not_ghost);
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
    facet_type[nb_facet_types++] = mesh.getFacetType(*it);
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

/* -------------------------------------------------------------------------- */
void MeshUtils::insertIntrinsicCohesiveElements(Mesh & mesh,
						Mesh & mesh_facets,
						ElementType type_facet,
						const Array<bool> & facet_insertion) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  /// create dummy doubled data containers
  ByElementTypeUInt doubled_facets("doubled_facets", "");
  mesh_facets.initByElementTypeArray(doubled_facets, 2, spatial_dimension - 1);

  /// setup byelementtype insertion data
  ByElementTypeArray<bool> facet_ins("facet_ins", "");
  mesh_facets.initByElementTypeArray(facet_ins, 1, spatial_dimension - 1);

  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1);

  for(; it != end; ++it) {

    if (*it != type_facet) continue;

    Array<bool> & f_ins = facet_ins(type_facet);
    f_ins.copy(facet_insertion);
  }

  /// add cohesive connectivity type
  ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
  mesh.addConnectivityType(type_cohesive);

  /// insert cohesive elements
  insertCohesiveElements(mesh,
			 mesh_facets,
			 facet_ins,
			 doubled_facets,
			 false);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::insertCohesiveElements(Mesh & mesh,
				       Mesh & mesh_facets,
				       const ByElementTypeArray<bool> & facet_insertion,
				       ByElementTypeUInt & doubled_facets,
				       const bool extrinsic) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  NewNodesEvent node_event;
  Array<UInt> & doubled_nodes = node_event.getList();
  doubled_nodes.extendComponentsInterlaced(2, 1);

  NewElementsEvent element_event;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      Array<UInt> & doubled_f = doubled_facets(type_facet, gt_facet);
      doubled_f.resize(0);

      const Array<bool> & f_insertion = facet_insertion(type_facet, gt_facet);
      Array<UInt> facets_to_be_doubled;

      for (UInt f = 0; f < f_insertion.getSize(); ++f) {
	if (f_insertion(f))
	  facets_to_be_doubled.push_back(f);
      }

      if(facets_to_be_doubled.getSize() == 0) continue;

      Element facet(type_facet, 0, gt_facet);

      /// update mesh
      for (UInt f = 0; f < facets_to_be_doubled.getSize(); ++f) {
	facet.element = facets_to_be_doubled(f);
	doubleFacet(mesh, mesh_facets, facet, doubled_nodes, doubled_f);
      }

      /// double middle nodes if it's the case
      if (type_facet == _segment_3)
	doubleMiddleNode(mesh, mesh_facets, gt_facet, doubled_nodes, doubled_f);

      /// loop over doubled facets to insert cohesive elements
      ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
      Array<UInt> & conn_cohesive = mesh.getConnectivity(type_cohesive, gt_facet);
      const Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet,
								   gt_facet);

      UInt nb_nodes_per_facet = conn_facet.getNbComponent();
      const Array<Real> & position = mesh.getNodes();
      Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(type_facet, gt_facet);

      const Real epsilon = std::numeric_limits<Real>::epsilon();

      UInt old_nb_cohesive_elements = conn_cohesive.getSize();
      conn_cohesive.resize(old_nb_cohesive_elements + doubled_f.getSize());

      for (UInt f = 0; f < doubled_f.getSize(); ++f) {
	UInt nb_cohesive_elements = old_nb_cohesive_elements + f;

	UInt first_facet = doubled_f(f, 0);
	UInt second_facet = doubled_f(f, 1);

	/// copy first facet's connectivity
	for (UInt n = 0; n < nb_nodes_per_facet; ++n)
	  conn_cohesive(nb_cohesive_elements, n) = conn_facet(first_facet, n);

	/// check if first nodes of the two facets are coincident or not
	UInt first_facet_node = conn_facet(first_facet, 0);
	UInt second_facet_node = conn_facet(second_facet, 0);
	bool second_facet_inversion = false;

	for (UInt dim = 0; dim < mesh.getSpatialDimension(); ++dim) {
	  if (std::abs( (position(first_facet_node, dim)
			 - position(second_facet_node, dim))
			/ position(second_facet_node, dim)) >= epsilon) {
	    second_facet_inversion = true;
	    break;
	  }
	}

	/// if the two nodes are coincident second facet connectivity is
	/// normally copied, otherwise the copy is reverted
	if (!second_facet_inversion) {
	  for (UInt n = 0; n < nb_nodes_per_facet; ++n)
	    conn_cohesive(nb_cohesive_elements, n + nb_nodes_per_facet)
	      = conn_facet(second_facet, n);
	}
	else {
	  for (UInt n = 0; n < nb_nodes_per_facet; ++n)
	    conn_cohesive(nb_cohesive_elements, n + nb_nodes_per_facet)
	      = conn_facet(second_facet, nb_nodes_per_facet - n - 1);
	}

	/// update element_to_facet vectors
	Element cohesive_element(type_cohesive, nb_cohesive_elements,
				 gt_facet, _ek_cohesive);
	element_to_facet(first_facet)[1] = cohesive_element;
	element_to_facet(second_facet)[1] = cohesive_element;
      }
    }
  }

  UInt total_nb_new_nodes = 0;
  UInt nb_new_nodes = doubled_nodes.getSize();

  /// update node global id
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if (psize == 1) total_nb_new_nodes = nb_new_nodes;
  else {
    if (extrinsic) {
      Array<Int> & nodes_type = const_cast<Array<Int> &>(mesh.getNodesType());
      UInt nb_old_nodes = nodes_type.getSize();
      nodes_type.resize(nb_old_nodes + nb_new_nodes);

      /// update nodes' type
      for (UInt n = 0; n < nb_new_nodes; ++n) {
	UInt old_node = doubled_nodes(n, 0);
	UInt new_node = doubled_nodes(n, 1);
	nodes_type(new_node) = nodes_type(old_node);
      }

      Int * nb_nodes_per_proc = new Int[psize];
      nb_nodes_per_proc[prank] = nb_new_nodes;
      comm.allGather(nb_nodes_per_proc, 1);

      Array<UInt> & nodes_global_ids =
	const_cast<Array<UInt> &>(mesh.getGlobalNodesIds());

      nodes_global_ids.resize(nb_old_nodes + nb_new_nodes);
      UInt nb_global_nodes = mesh.getNbGlobalNodes();
      UInt current_added_node = nb_global_nodes;

      for (Int p = 0; p < prank; ++p) {
	current_added_node += nb_nodes_per_proc[p];
	total_nb_new_nodes += nb_nodes_per_proc[p];
      }

      for (Int p = prank; p < psize; ++p)
	total_nb_new_nodes += nb_nodes_per_proc[p];

      for (UInt n = 0; n < nb_new_nodes; ++n, ++current_added_node)
	nodes_global_ids(nb_old_nodes + n) = current_added_node;
    }
  }

  if (total_nb_new_nodes > 0 || !extrinsic) {
    mesh.updateTypesOffsets(_not_ghost);
    mesh.sendEvent(node_event);
    mesh.sendEvent(element_event);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doubleMiddleNode(Mesh & mesh,
				 Mesh & mesh_facets,
				 GhostType gt_facet,
				 Array<UInt> & doubled_nodes,
				 const Array<UInt> & doubled_facets) {

  AKANTU_DEBUG_IN();

  ElementType type_facet = _segment_3;

  Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
  Array<Real> & position = mesh.getNodes();
  UInt spatial_dimension = mesh.getSpatialDimension();

  const Array< std::vector<Element> > & elem_to_facet
    = mesh_facets.getElementToSubelement(type_facet, gt_facet);

  UInt nb_new_facets = doubled_facets.getSize();

  UInt old_nb_doubled_nodes = doubled_nodes.getSize();
  doubled_nodes.resize(old_nb_doubled_nodes + nb_new_facets);

  UInt old_nb_nodes = position.getSize();
  position.resize(old_nb_nodes + nb_new_facets);

  for (UInt f = 0; f < nb_new_facets; ++f) {

    UInt facet_first = doubled_facets(f, 0);
    UInt facet_second = doubled_facets(f, 1);

    UInt new_node = old_nb_nodes + f;

    /// store doubled nodes
    UInt nb_doubled_nodes = old_nb_doubled_nodes + f;

    UInt old_node = conn_facet(facet_first, 2);

    doubled_nodes(nb_doubled_nodes, 0) = old_node;
    doubled_nodes(nb_doubled_nodes, 1) = new_node;

    /// update position
    for (UInt dim = 0; dim < spatial_dimension; ++dim)
      position(new_node, dim) = position(old_node, dim);

    /// update facet connectivity
    conn_facet(facet_second, 2) = new_node;

    /// update element connectivity
    for (UInt el = 0; el < elem_to_facet(facet_second).size(); ++el) {
      Element elem = elem_to_facet(facet_second)[el];
      ElementType type_elem = elem.type;
      if (type_elem != _not_defined) {
	UInt elem_global = elem.element;
	GhostType gt_elem = elem.ghost_type;
	Array<UInt> & conn_elem = mesh.getConnectivity(type_elem, gt_elem);

	for (UInt n = 0; n < conn_elem.getNbComponent(); ++n) {
	  if (conn_elem(elem_global, n) == old_node)
	    conn_elem(elem_global, n) = new_node;
	}
      }
    }

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doubleFacet(Mesh & mesh,
			    Mesh & mesh_facets,
			    Element & facet,
			    Array<UInt> & doubled_nodes,
			    Array<UInt> & doubled_facets) {
  AKANTU_DEBUG_IN();

  const UInt f_index = facet.element;
  const ElementType type_facet = facet.type;
  const GhostType gt_facet = facet.ghost_type;

  const ElementType type_subfacet = mesh.getFacetType(type_facet);
  const UInt nb_subfacet = mesh.getNbFacetsPerElement(type_facet);

  Array< std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet, gt_facet);

  if (element_to_facet(f_index)[1].type == _not_defined ||
      element_to_facet(f_index)[1].kind == _ek_cohesive) {
    AKANTU_DEBUG_WARNING("attempt to double a facet on the boundary");
    return;
  }

  Array<Element> & subfacet_to_facet
    = mesh_facets.getSubelementToElement(type_facet, gt_facet);

  /// adding a new facet by copying original one

  /// create new connectivity
  Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
  UInt nb_facet = conn_facet.getSize();
  conn_facet.resize(nb_facet + 1);
  for (UInt n = 0; n < conn_facet.getNbComponent(); ++n)
    conn_facet(nb_facet, n) = conn_facet(f_index, n);

  /// store doubled facets
  UInt nb_doubled_facets = doubled_facets.getSize();
  doubled_facets.resize(nb_doubled_facets + 1);
  doubled_facets(nb_doubled_facets, 0) = f_index;
  doubled_facets(nb_doubled_facets, 1) = nb_facet;

  /// update elements connected to facet
  std::vector<Element> first_facet_list = element_to_facet(f_index);
  element_to_facet.push_back(first_facet_list);

  /// set new and original facets as boundary facets
  element_to_facet(f_index)[1] = ElementNull;

  element_to_facet(nb_facet)[0] = element_to_facet(nb_facet)[1];
  element_to_facet(nb_facet)[1] = ElementNull;

  /// update facet_to_element vector
  ElementType type = element_to_facet(nb_facet)[0].type;
  UInt el = element_to_facet(nb_facet)[0].element;
  GhostType gt = element_to_facet(nb_facet)[0].ghost_type;
  Array<Element> & facet_to_element = mesh_facets.getSubelementToElement(type, gt);

  UInt i;
  for (i = 0; facet_to_element(el, i) != facet
	 && i <= facet_to_element.getNbComponent(); ++i);

  facet_to_element(el, i).element = nb_facet;

  /// create new links to subfacets and update list of facets
  /// connected to subfacets
  subfacet_to_facet.resize(nb_facet + 1);
  for (UInt sf = 0; sf < nb_subfacet; ++sf) {
    subfacet_to_facet(nb_facet, sf) = subfacet_to_facet(f_index, sf);

    Element subfacet = subfacet_to_facet(f_index, sf);
    if (subfacet == ElementNull) continue;
    UInt sf_index = subfacet.element;
    GhostType sf_gt = subfacet.ghost_type;

    Array< std::vector<Element> > & facet_to_subfacet
      = mesh_facets.getElementToSubelement(type_subfacet, sf_gt);

    /// find index to start looping around facets connected to current
    /// subfacet
    UInt start = 0;
    UInt nb_connected_facets = facet_to_subfacet(sf_index).size();

    while (facet_to_subfacet(sf_index)[start] != facet
    	   && start <= facet_to_subfacet(sf_index).size()) ++start;

    /// add the new facet to the list next to the original one
    ++nb_connected_facets;
    facet_to_subfacet(sf_index).resize(nb_connected_facets);

    for (UInt f = nb_connected_facets - 1; f > start; --f) {
      facet_to_subfacet(sf_index)[f] = facet_to_subfacet(sf_index)[f - 1];
    }


    /// check if the new facet should be inserted before or after the
    /// original one: the second element connected to both original
    /// and new facet will be _not_defined, so I check if the first
    /// one is equal to one of the elements connected to the following
    /// facet in the facet_to_subfacet vector
    UInt f_start = facet_to_subfacet(sf_index)[start].element;
    Element facet_next;
    if (start + 2 == nb_connected_facets)
      facet_next = facet_to_subfacet(sf_index)[0];
    else
      facet_next = facet_to_subfacet(sf_index)[start + 2];

    Array< std::vector<Element> > & element_to_facet_next
      = mesh_facets.getElementToSubelement(facet_next.type, facet_next.ghost_type);

    UInt f_next = facet_next.element;

    if ((element_to_facet(f_start)[0] == element_to_facet_next(f_next)[0])
	|| ( element_to_facet(f_start)[0] == element_to_facet_next(f_next)[1]))
      facet_to_subfacet(sf_index)[start].element = nb_facet;
    else
      facet_to_subfacet(sf_index)[start + 1].element = nb_facet;

    /// loop on every facet connected to the current subfacet
    for (UInt f = start + 2; ; ++f) {

      /// reset f in order to continue looping from the beginning
      if (f == nb_connected_facets) f = 0;
      /// exit loop if it reaches the end
      if (f == start) break;

      /// if current loop facet is on the boundary, double subfacet
      Element f_global = facet_to_subfacet(sf_index)[f];

      Array< std::vector<Element> > & el_to_f
	= mesh_facets.getElementToSubelement(f_global.type, f_global.ghost_type);

      if (el_to_f(f_global.element)[1].type == _not_defined ||
	  el_to_f(f_global.element)[1].kind == _ek_cohesive) {
	doubleSubfacet(mesh,
		       mesh_facets,
		       subfacet,
		       start,
		       f,
		       doubled_nodes);
	break;
      }

    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doubleSubfacet(Mesh & mesh,
			       Mesh & mesh_facets,
			       const Element & subfacet,
			       UInt start,
			       UInt end,
			       Array<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  const UInt sf_index = subfacet.element;
  const ElementType type_subfacet = subfacet.type;
  const GhostType gt_subfacet = subfacet.ghost_type;

  Array< std::vector<Element> > & facet_to_subfacet
    = mesh_facets.getElementToSubelement(type_subfacet, gt_subfacet);
  UInt nb_subfacet = facet_to_subfacet.getSize();

  Array<UInt> & conn_point = mesh_facets.getConnectivity(_point_1, gt_subfacet);
  Array<Real> & position = mesh.getNodes();
  UInt spatial_dimension = mesh.getSpatialDimension();

  /// add the new subfacet
  if (spatial_dimension == 2) {

    UInt new_node = position.getSize();

    /// add new node in connectivity
    UInt new_subfacet = conn_point.getSize();
    conn_point.resize(new_subfacet + 1);
    conn_point(new_subfacet) = new_node;

    /// store doubled nodes
    UInt nb_doubled_nodes = doubled_nodes.getSize();
    doubled_nodes.resize(nb_doubled_nodes + 1);

    UInt old_node = doubled_nodes(nb_doubled_nodes, 0)
      = conn_point(sf_index);
    doubled_nodes(nb_doubled_nodes, 1) = new_node;

    /// update position
    position.resize(new_node + 1);
    for (UInt dim = 0; dim < spatial_dimension; ++dim)
      position(new_node, dim) = position(old_node, dim);
  }

  /// create a vector for the new subfacet in facet_to_subfacet
  facet_to_subfacet.resize(nb_subfacet + 1);

  UInt nb_connected_facets = facet_to_subfacet(sf_index).size();
  UInt new_node = conn_point(conn_point.getSize() - 1);
  UInt old_node = conn_point(sf_index);
  /// loop over facets from start to end
  for (UInt f = start + 1; ; ++f) {

    /// reset f in order to continue looping from the beginning
    if (f == nb_connected_facets) f = 0;

    Element facet_el = facet_to_subfacet(sf_index)[f];
    UInt f_global = facet_el.element;
    ElementType type_facet = facet_el.type;
    GhostType gt_facet = facet_el.ghost_type;

    Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
    UInt nb_nodes_per_facet = conn_facet.getNbComponent();

    /// update facet connectivity
    UInt i;
    for (i = 0; conn_facet(f_global, i) != old_node
	   && i <= nb_nodes_per_facet; ++i);
    conn_facet(f_global, i) = new_node;
    UInt facet_conn_node = i;

    /// update element connectivity
    const Array< std::vector<Element> > & elem_to_facet
      = mesh_facets.getElementToSubelement(type_facet, gt_facet);

    for (UInt el = 0; el < elem_to_facet(f_global).size(); ++el) {
      Element elem = elem_to_facet(f_global)[el];
      ElementType type_elem = elem.type;
      if (type_elem != _not_defined) {
	UInt elem_global = elem.element;
	GhostType gt_elem = elem.ghost_type;
	Array<UInt> & conn_elem = mesh.getConnectivity(type_elem, gt_elem);
	UInt nb_nodes_per_element = conn_elem.getNbComponent();

	/// integer to do the final for loop
	UInt n = 0;

	/// cohesive elements: it's neccesary to identify the correct
	/// facet to be updated by looking for another node
	if (elem_to_facet(f_global)[el].kind == _ek_cohesive) {

	  UInt second_node = 0;
	  if (facet_conn_node == 0) second_node = 1;

	  UInt conn_second_node = conn_facet(f_global, second_node);
	  for (i = 0; i < nb_nodes_per_element &&
		 conn_elem(elem_global, i) != conn_second_node; ++i);

	  if (i >= nb_nodes_per_facet) n = nb_nodes_per_facet;
	}
	/// update connectivity
	for (; n < nb_nodes_per_element; ++n) {
	  if (conn_elem(elem_global, n) == old_node) {
	    conn_elem(elem_global, n) = new_node;
	    break;
	  }
	}
      }
    }

    /// update subfacet_to_facet vector
    Array<Element> & subfacet_to_facet
      = mesh_facets.getSubelementToElement(type_facet, gt_facet);

    for (i = 0; subfacet_to_facet(f_global, i) != subfacet
	   && i <= subfacet_to_facet.getNbComponent(); ++i);
    subfacet_to_facet(f_global, i).element = nb_subfacet;

    /// add current facet to facet_to_subfacet last position
    Element current_facet(type_facet, f_global, gt_facet);
    facet_to_subfacet(nb_subfacet).push_back(current_facet);

    /// exit loop if it reaches the end
    if (f == end) break;
  }

  /// rearrange the facets connected to the original subfacet and
  /// compute the new number of facets connected to it
  if (end < start) {
    for (UInt f = 0; f < start - end; ++f)
      facet_to_subfacet(sf_index)[f] = facet_to_subfacet(sf_index)[f + end + 1];

    nb_connected_facets = start - end;
  }
  else {
    for (UInt f = 1; f < nb_connected_facets - end; ++f)
      facet_to_subfacet(sf_index)[start + f] = facet_to_subfacet(sf_index)[end + f];

    nb_connected_facets -= end - start;
  }

  /// resize list of facets of the original subfacet
  facet_to_subfacet(sf_index).resize(nb_connected_facets);

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__

//  LocalWords:  ElementType
