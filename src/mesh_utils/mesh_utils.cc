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
  if (spatial_dimension == _all_dimensions) spatial_dimension = mesh.getSpatialDimension();


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
  if (spatial_dimension == _all_dimensions) spatial_dimension = mesh.getSpatialDimension();
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
						CSR<UInt> & node_to_elem,
						const ElementType & type,
						const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  UInt nb_nodes = mesh.getNbNodes();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_elements = mesh.getConnectivity(type, ghost_type).getSize();

  UInt * conn_val = mesh.getConnectivity(type, ghost_type).storage();

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
  ByElementTypeReal barycenter("barycenter_tmp", mesh.getID());
  ByElementTypeUInt prank_to_element("prank_to_elem", mesh.getID());

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType gt_facet = *gt;
    Mesh::type_iterator it  = mesh.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension - 1, gt_facet);
    for(; it != end; ++it) {
      mesh.getConnectivity(*it, *gt).resize(0);
      // \todo inform the mesh event handler
    }
  }

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

  ByElementTypeUInt prank_to_element("prank_to_elem", mesh.getID());
  buildAllFacetsParallel(mesh, mesh_facets, prank_to_element);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::buildAllFacetsParallel(Mesh & mesh,
				       Mesh & mesh_facets,
				       ByElementTypeUInt & prank_to_element) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  ByElementTypeReal barycenter("barycenter_tmp", mesh.getID());

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

  /// copy nodes type pointer
  mesh_facets.nodes_type = mesh.nodes_type;

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
      Array<UInt> * connectivity_facets =
	mesh_facets.getConnectivityPointer(facet_type, ghost_type);

      connectivity_facets->resize(0);

      Array< std::vector<Element> > * element_to_subelement =
	mesh_facets.getElementToSubelementPointer(facet_type, ghost_type);

      element_to_subelement->resize(0);

      Array<Element> * subelement_to_element =
	mesh_facets.getSubelementToElementPointer(type, ghost_type);

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
      Array< std::vector<Element> > * element_to_subelement =
	&mesh_facets.getElementToSubelement(facet_type, ghost_type);
      Array<UInt> * connectivity_facets =
	&mesh_facets.getConnectivity(facet_type, ghost_type);

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
	  for (; first_node_elements != first_node_elements_end;
	       ++first_node_elements, ++local_el) {
	    for (UInt n = 1; n < nb_nodes_per_facet; ++n) {
	      CSR<Element>::iterator node_elements_begin = node_to_elem.begin(facet(n));
	      CSR<Element>::iterator node_elements_end   = node_to_elem.end  (facet(n));
	      counter(local_el) += std::count(node_elements_begin,
					      node_elements_end,
					      *first_node_elements);
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
		  if (dimension == spatial_dimension - 1)
                    MeshUtils::sortElements(elements, facet, mesh, mesh_facets, barycenter);
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

    Mesh::type_iterator it  = mesh.firstType(_all_dimensions, ghost_type, _ek_not_defined);
    Mesh::type_iterator end = mesh.lastType(_all_dimensions, ghost_type, _ek_not_defined);
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
void MeshUtils::insertIntrinsicCohesiveElementsInArea(Mesh & mesh,
						      const Array<Real> & limits) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(limits.getNbComponent() == 2,
		      "Number of components for limits array must be 2");

  UInt spatial_dimension = mesh.getSpatialDimension();

  AKANTU_DEBUG_ASSERT(limits.getSize() == spatial_dimension,
		      "Limits array size must be equal to spatial dimension");

  Real * bary_facet = new Real[spatial_dimension];

  Mesh mesh_facets(spatial_dimension, mesh.getNodes(), "mesh_facets");
  buildAllFacets(mesh, mesh_facets);

  resetFacetToDouble(mesh_facets);

  ByElementTypeArray<bool> facet_insertion("facet_insertion", "");
  for (UInt sp = 0; sp < spatial_dimension - 1; ++sp)
    mesh_facets.initByElementTypeArray(facet_insertion, 1, sp - 1,
				       false, _ek_regular, true);

  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1);

  for(; it != end; ++it) {
    const ElementType type_facet = *it;
    Array<bool> & f_insertion = facet_insertion(type_facet);
    Array<std::vector<Element> > & element_to_facet
      = mesh_facets.getElementToSubelement(type_facet);

    UInt nb_facet = mesh_facets.getNbElement(type_facet);

    for (UInt f = 0; f < nb_facet; ++f) {

      if (element_to_facet(f)[1] == ElementNull) continue;

      mesh_facets.getBarycenter(f, type_facet, bary_facet);

      UInt coord_in_limit = 0;

      while (coord_in_limit < spatial_dimension &&
	     bary_facet[coord_in_limit] > limits(coord_in_limit, 0) &&
	     bary_facet[coord_in_limit] < limits(coord_in_limit, 1))
	++coord_in_limit;

      if (coord_in_limit == spatial_dimension)
	f_insertion(f) = true;
    }
  }

  delete [] bary_facet;

  insertCohesiveElements(mesh,
			 mesh_facets,
			 facet_insertion,
			 false);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::insertIntrinsicCohesiveElements(Mesh & mesh,
						Mesh & mesh_facets,
						ElementType type_facet,
						const Array<bool> & facet_insertion) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  resetFacetToDouble(mesh_facets);

  /// add cohesive connectivity type
  ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
  mesh.addConnectivityType(type_cohesive);

  /// setup byelementtype insertion data
  ByElementTypeArray<bool> facet_ins("facet_ins", "");
  for (UInt sp = 0; sp < spatial_dimension - 1; ++sp)
    mesh_facets.initByElementTypeArray(facet_ins, 1, sp - 1,
				       false, _ek_regular, true);

  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1);

  for(; it != end; ++it) {

    if (*it != type_facet) continue;

    Array<bool> & f_ins = facet_ins(type_facet);
    f_ins.copy(facet_insertion);
  }

  /// insert cohesive elements
  insertCohesiveElements(mesh,
			 mesh_facets,
			 facet_ins,
			 false);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::insertCohesiveElements(Mesh & mesh,
				       Mesh & mesh_facets,
				       ByElementTypeArray<bool> & facet_insertion,
				       const bool extrinsic) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  NewNodesEvent node_event;
  Array<UInt> & doubled_nodes = node_event.getList();
  doubled_nodes.extendComponentsInterlaced(2, 1);

  UInt total_nb_new_elements = 0;
  NewElementsEvent element_event;

  if (updateFacetToDouble(mesh_facets, facet_insertion)) {

    if (spatial_dimension == 1) {
      doublePointFacet(mesh, mesh_facets, doubled_nodes);
    }
    else {
      doubleFacet(mesh, mesh_facets, spatial_dimension - 1,
		  doubled_nodes, true);
      findSubfacetToDouble<false>(mesh, mesh_facets);

      if (spatial_dimension == 2) {
	doubleSubfacet<2>(mesh, mesh_facets, doubled_nodes);
      }
      else if (spatial_dimension == 3) {
	doubleFacet(mesh, mesh_facets, 1, doubled_nodes, false);
	findSubfacetToDouble<true>(mesh, mesh_facets);
	doubleSubfacet<3>(mesh, mesh_facets, doubled_nodes);
      }
    }

    updateCohesiveData(mesh, mesh_facets, element_event);

    for (ghost_type_t::iterator gt = ghost_type_t::begin();
	 gt != ghost_type_t::end(); ++gt) {

      GhostType gt_facet = *gt;

      Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1,
						      gt_facet);
      Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1,
						     gt_facet);

      for(; it != end; ++it) {
	ElementType type_facet = *it;

	Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							      type_facet,
							      gt_facet);
	total_nb_new_elements += f_to_double.getSize();

	/// reset facet insertion
	facet_insertion(type_facet, gt_facet).clear();
      }
    }

    resetFacetToDouble(mesh_facets);
  }

  UInt nb_new_nodes = doubled_nodes.getSize();
  UInt total_nb_new_nodes = nb_new_nodes;

  /// @todo update node global id

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();

  if (psize > 1 && extrinsic) {
    Array<Int> & nodes_type = *mesh.nodes_type;
    UInt nb_old_nodes = nodes_type.getSize();
    nodes_type.resize(nb_old_nodes + nb_new_nodes);

    /// update nodes' type
    for (UInt n = 0; n < nb_new_nodes; ++n) {
      UInt old_node = doubled_nodes(n, 0);
      UInt new_node = doubled_nodes(n, 1);
      nodes_type(new_node) = nodes_type(old_node);
    }

    comm.allReduce(&total_nb_new_nodes, 1, _so_sum);
    comm.allReduce(&total_nb_new_elements, 1, _so_sum);

    (*mesh.nodes_global_ids).resize(nb_old_nodes + nb_new_nodes);
  }

  mesh.nb_global_nodes += total_nb_new_nodes;
  mesh_facets.nb_global_nodes += total_nb_new_nodes;

  if (total_nb_new_nodes > 0 || !extrinsic) {
    mesh.sendEvent(node_event);
  }

  if (total_nb_new_elements > 0 || !extrinsic) {
    mesh.updateTypesOffsets(_not_ghost);
    mesh.sendEvent(element_event);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doubleNodes(Mesh & mesh,
			    const Array<UInt> & old_nodes,
			    Array<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  Array<Real> & position = mesh.getNodes();
  UInt spatial_dimension = mesh.getSpatialDimension();

  UInt old_nb_nodes = position.getSize();
  UInt new_nb_nodes = old_nb_nodes + old_nodes.getSize();

  UInt old_nb_doubled_nodes = doubled_nodes.getSize();
  UInt new_nb_doubled_nodes = old_nb_doubled_nodes + old_nodes.getSize();

  position.resize(new_nb_nodes);
  doubled_nodes.resize(new_nb_doubled_nodes);

  for (UInt n = 0; n < old_nodes.getSize(); ++n) {
    UInt new_node = old_nb_nodes + n;

    /// store doubled nodes
    doubled_nodes(old_nb_doubled_nodes + n, 0) = old_nodes(n);
    doubled_nodes(old_nb_doubled_nodes + n, 1) = new_node;

    /// update position
    std::copy(position.begin(spatial_dimension) + old_nodes(n),
	      position.begin(spatial_dimension) + old_nodes(n) + 1,
	      position.begin(spatial_dimension) + new_node);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doubleFacet(Mesh & mesh,
			    Mesh & mesh_facets,
			    UInt facet_dimension,
			    Array<UInt> & doubled_nodes,
			    bool facet_mode) {
  AKANTU_DEBUG_IN();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(facet_dimension, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(facet_dimension, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							    type_facet,
							    gt_facet);
      UInt nb_facet_to_double = f_to_double.getSize();

      if (nb_facet_to_double == 0) continue;

      ElementType type_subfacet = Mesh::getFacetType(type_facet);
      const UInt nb_subfacet_per_facet = Mesh::getNbFacetsPerElement(type_facet);
      GhostType gt_subfacet = _casper;
      Array<std::vector<Element> > * f_to_subfacet = NULL;

      Array<Element> & subfacet_to_facet
	= mesh_facets.getSubelementToElement(type_facet, gt_facet);

      Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
      UInt nb_nodes_per_facet = conn_facet.getNbComponent();

      UInt old_nb_facet = conn_facet.getSize();
      UInt new_nb_facet = old_nb_facet + nb_facet_to_double;

      conn_facet.resize(new_nb_facet);
      subfacet_to_facet.resize(new_nb_facet);

      UInt new_facet = old_nb_facet - 1;
      Element new_facet_el(type_facet, 0, gt_facet);

      for (UInt facet = 0; facet < nb_facet_to_double; ++facet) {
	UInt old_facet = f_to_double(facet);
	++new_facet;

	/// adding a new facet by copying original one

	/// copy connectivity in new facet
	std::copy(conn_facet.begin(nb_nodes_per_facet) + old_facet,
		  conn_facet.begin(nb_nodes_per_facet) + old_facet + 1,
		  conn_facet.begin(nb_nodes_per_facet) + new_facet);

	/// update subfacet_to_facet
	std::copy(subfacet_to_facet.begin(nb_subfacet_per_facet) + old_facet,
		  subfacet_to_facet.begin(nb_subfacet_per_facet) + old_facet + 1,
		  subfacet_to_facet.begin(nb_subfacet_per_facet) + new_facet);

	new_facet_el.element = new_facet;

	/// loop on every subfacet
	for (UInt sf = 0; sf < nb_subfacet_per_facet; ++sf) {
	  Element & subfacet = subfacet_to_facet(old_facet, sf);
	  if (subfacet == ElementNull) continue;

	  if (gt_subfacet != subfacet.ghost_type) {
	    gt_subfacet = subfacet.ghost_type;
	    f_to_subfacet = & mesh_facets.getElementToSubelement(type_subfacet,
								 subfacet.ghost_type);
	  }

	  /// update facet_to_subfacet array
	  (*f_to_subfacet)(subfacet.element).push_back(new_facet_el);
	}
      }

      /// update facet_to_subfacet and _segment_3 facets if any
      if (!facet_mode) {
	updateSubfacetToFacet(mesh_facets, type_facet, gt_facet, true);
	updateFacetToSubfacet(mesh_facets, type_facet, gt_facet, true);
	updateQuadraticSegments<true>(mesh, mesh_facets, type_facet,
				      gt_facet, doubled_nodes);
      }
      else
	updateQuadraticSegments<false>(mesh, mesh_facets, type_facet,
				       gt_facet, doubled_nodes);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool MeshUtils::updateFacetToDouble(Mesh & mesh_facets,
				    const ByElementTypeArray<bool> & facet_insertion) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();
  bool facets_to_double = false;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      const Array<bool> & f_insertion = facet_insertion(type_facet, gt_facet);
      Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							    type_facet,
							    gt_facet);

      Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(type_facet, gt_facet);

      ElementType el_type = _not_defined;
      GhostType el_gt = _casper;
      UInt nb_facet_per_element = 0;
      Element old_facet_el(type_facet, 0, gt_facet);

      Array<Element> * facet_to_element = NULL;

      for (UInt f = 0; f < f_insertion.getSize(); ++f) {

	if (f_insertion(f) == false) continue;

	facets_to_double = true;

	if (element_to_facet(f)[1].type == _not_defined ||
	    element_to_facet(f)[1].kind == _ek_cohesive) {
	  AKANTU_DEBUG_WARNING("attempt to double a facet on the boundary");
	  continue;
	}

	f_to_double.push_back(f);

	UInt new_facet =
	  mesh_facets.getNbElement(type_facet, gt_facet) + f_to_double.getSize() - 1;
	old_facet_el.element = f;

	/// update facet_to_element vector
	Element & elem_to_update = element_to_facet(f)[1];
	UInt el = elem_to_update.element;

	if (elem_to_update.ghost_type != el_gt || elem_to_update.type != el_type) {
	  el_type = elem_to_update.type;
	  el_gt = elem_to_update.ghost_type;
	  facet_to_element = & mesh_facets.getSubelementToElement(el_type, el_gt);
	  nb_facet_per_element = facet_to_element->getNbComponent();
	}

	Element * f_update
	  = std::find(facet_to_element->storage() + el * nb_facet_per_element,
		      facet_to_element->storage() + el * nb_facet_per_element
		      + nb_facet_per_element,
		      old_facet_el);

	AKANTU_DEBUG_ASSERT(facet_to_element->storage() + el * nb_facet_per_element !=
			    facet_to_element->storage() + el * nb_facet_per_element
			    + nb_facet_per_element,
			    "Facet not found");

	f_update->element = new_facet;

	/// update elements connected to facet
	std::vector<Element> first_facet_list = element_to_facet(f);
	element_to_facet.push_back(first_facet_list);

	/// set new and original facets as boundary facets
	element_to_facet(new_facet)[0] = element_to_facet(new_facet)[1];

	element_to_facet(f)[1] = ElementNull;
	element_to_facet(new_facet)[1] = ElementNull;
	
      }
    }
  }

  AKANTU_DEBUG_OUT();
  return facets_to_double;
}

/* -------------------------------------------------------------------------- */
void MeshUtils::resetFacetToDouble(Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  for(UInt g = _not_ghost; g <= _ghost; ++g) {
    GhostType gt = (GhostType) g;

    Mesh::type_iterator it  = mesh_facets.firstType(_all_dimensions, gt);
    Mesh::type_iterator end = mesh_facets.lastType(_all_dimensions, gt);
    for(; it != end; ++it) {
      ElementType type = *it;
      mesh_facets.getDataPointer<UInt>("facet_to_double", type, gt, 1)->resize(0);

      mesh_facets.getDataPointer<std::vector<Element> >
	("facets_to_subfacet_double", type, gt, 1)->resize(0);

      mesh_facets.getDataPointer<std::vector<Element> >
	("elements_to_subfacet_double", type, gt, 1)->resize(0);

      mesh_facets.getDataPointer<std::vector<Element> >
	("subfacets_to_subsubfacet_double", type, gt, 1)->resize(0);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool subsubfacet_mode>
void MeshUtils::findSubfacetToDouble(Mesh & mesh, Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();
  if (spatial_dimension == 1) return;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							    type_facet,
							    gt_facet);
      UInt nb_facet_to_double = f_to_double.getSize();
      if (nb_facet_to_double == 0) continue;

      ElementType type_subfacet = Mesh::getFacetType(type_facet);
      GhostType gt_subfacet = _casper;

      ElementType type_subsubfacet = Mesh::getFacetType(type_subfacet);
      GhostType gt_subsubfacet = _casper;

      Array<UInt> * conn_subfacet = NULL;
      Array<UInt> * sf_to_double = NULL;
      Array<std::vector<Element> > * sf_to_subfacet_double = NULL;
      Array<std::vector<Element> > * f_to_subfacet_double = NULL;
      Array<std::vector<Element> > * el_to_subfacet_double = NULL;

      UInt nb_subfacet = Mesh::getNbFacetsPerElement(type_facet);

      UInt nb_subsubfacet;
      UInt nb_nodes_per_sf_el;

      if (subsubfacet_mode) {
	nb_nodes_per_sf_el = Mesh::getNbNodesPerElement(type_subsubfacet);
	nb_subsubfacet = Mesh::getNbFacetsPerElement(type_subfacet);
      }
      else
	nb_nodes_per_sf_el = Mesh::getNbNodesPerElement(type_subfacet);

      Array<Element> & subfacet_to_facet
	= mesh_facets.getSubelementToElement(type_facet, gt_facet);

      Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(type_facet, gt_facet);

      Array<Element> * subsubfacet_to_subfacet = NULL;

      UInt old_nb_facet = subfacet_to_facet.getSize() - nb_facet_to_double;

      Element current_facet(type_facet, 0, gt_facet);
      std::vector<Element> element_list;
      std::vector<Element> facet_list;
      std::vector<Element> * subfacet_list;
      if (subsubfacet_mode)
	subfacet_list = new std::vector<Element>;

      /// map to filter subfacets
      Array< std::vector<Element> > * facet_to_subfacet = NULL;

      /// this is used to make sure that both new and old facets are
      /// checked
      UInt facets[2];

      /// loop on every facet
      for (UInt f_index = 0; f_index < 2; ++f_index) {
	for (UInt facet = 0; facet < nb_facet_to_double; ++facet) {
	  facets[bool(f_index)] = f_to_double(facet);
	  facets[!bool(f_index)] = old_nb_facet + facet;

	  UInt old_facet = facets[0];
	  UInt new_facet = facets[1];

	  Element & starting_element = element_to_facet(new_facet)[0];
	  current_facet.element = old_facet;

	  /// loop on every subfacet
	  for (UInt sf = 0; sf < nb_subfacet; ++sf) {

	    Element & subfacet = subfacet_to_facet(old_facet, sf);
	    if (subfacet == ElementNull) continue;

	    if (gt_subfacet != subfacet.ghost_type) {
	      gt_subfacet = subfacet.ghost_type;

	      if (subsubfacet_mode) {
		subsubfacet_to_subfacet
		  = & mesh_facets.getSubelementToElement(type_subfacet,
							 gt_subfacet);
	      }
	      else {
		conn_subfacet = & mesh_facets.getConnectivity(type_subfacet,
							      gt_subfacet);
		sf_to_double = & mesh_facets.getData<UInt>("facet_to_double",
							   type_subfacet,
							   gt_subfacet);

		f_to_subfacet_double
		  = & mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
								 type_subfacet,
								 gt_subfacet);

		el_to_subfacet_double
		  = & mesh_facets.getData<std::vector<Element> >("elements_to_subfacet_double",
								 type_subfacet,
								 gt_subfacet);

		facet_to_subfacet
		  = & mesh_facets.getElementToSubelement(type_subfacet,
							 gt_subfacet);
	      }
	    }

	    if (subsubfacet_mode) {
	      /// loop on every subsubfacet
	      for (UInt ssf = 0; ssf < nb_subsubfacet; ++ssf) {
		Element & subsubfacet
		  = (*subsubfacet_to_subfacet)(subfacet.element, ssf);

		if (subsubfacet == ElementNull) continue;

		if (gt_subsubfacet != subsubfacet.ghost_type) {
		  gt_subsubfacet = subsubfacet.ghost_type;
		  conn_subfacet = & mesh_facets.getConnectivity(type_subsubfacet,
								gt_subsubfacet);
		  sf_to_double = & mesh_facets.getData<UInt>("facet_to_double",
							     type_subsubfacet,
							     gt_subsubfacet);

		  sf_to_subfacet_double
		    = & mesh_facets.getData<std::vector<Element> >("subfacets_to_subsubfacet_double",
								   type_subsubfacet,
								   gt_subsubfacet);

		  f_to_subfacet_double
		    = & mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
								   type_subsubfacet,
								   gt_subsubfacet);

		  el_to_subfacet_double
		    = & mesh_facets.getData<std::vector<Element> >("elements_to_subfacet_double",
								   type_subsubfacet,
								   gt_subsubfacet);

		  facet_to_subfacet
		    = & mesh_facets.getElementToSubelement(type_subsubfacet,
							   gt_subsubfacet);
		}

		UInt global_ssf = subsubfacet.element;

		Vector<UInt> subsubfacet_connectivity(conn_subfacet->storage()
						      + global_ssf
						      * nb_nodes_per_sf_el,
						      nb_nodes_per_sf_el);

		/// check if subsubfacet is to be doubled
		if (findElementsAroundSubfacet<true>(mesh, mesh_facets,
						     starting_element, current_facet,
						     subsubfacet_connectivity,
						     element_list,
						     facet_list,
						     subfacet_list) == false &&
		    removeElementsInVector(*subfacet_list,
					   (*facet_to_subfacet)(global_ssf)) == false) {

		  sf_to_double->push_back(global_ssf);
		  sf_to_subfacet_double->push_back(*subfacet_list);
		  f_to_subfacet_double->push_back(facet_list);
		  el_to_subfacet_double->push_back(element_list);
		}
	      }
	    }
	    else {
	      const UInt global_sf = subfacet.element;

	      Vector<UInt> subfacet_connectivity(conn_subfacet->storage()
						 + global_sf
						 * nb_nodes_per_sf_el,
						 nb_nodes_per_sf_el);

	      /// check if subfacet is to be doubled
	      if (findElementsAroundSubfacet<false>(mesh, mesh_facets,
						    starting_element, current_facet,
						    subfacet_connectivity,
						    element_list,
						    facet_list) == false &&
		  removeElementsInVector(facet_list,
					 (*facet_to_subfacet)(global_sf)) == false) {

		sf_to_double->push_back(global_sf);
		f_to_subfacet_double->push_back(facet_list);
		el_to_subfacet_double->push_back(element_list);
	      }
	    }
	  }
	}
      }

      if (subsubfacet_mode)
	delete subfacet_list;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::updateCohesiveData(Mesh & mesh,
				   Mesh & mesh_facets,
				   NewElementsEvent & element_event) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  bool third_dimension = spatial_dimension == 3;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							    type_facet,
							    gt_facet);

      UInt nb_facet_to_double = f_to_double.getSize();
      if (nb_facet_to_double == 0) continue;

      ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
      Array<Element> * facet_to_coh_element =
	mesh.getSubelementToElementPointer(type_cohesive, gt_facet);

      Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
      Array<UInt> * conn_cohesive = mesh.getConnectivityPointer(type_cohesive,
								gt_facet);
      UInt nb_nodes_per_facet = Mesh::getNbNodesPerElement(type_facet);

      Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(type_facet, gt_facet);

      UInt old_nb_cohesive_elements = conn_cohesive->getSize();
      UInt new_nb_cohesive_elements = conn_cohesive->getSize() + nb_facet_to_double;

      UInt old_nb_facet = element_to_facet.getSize() - nb_facet_to_double;
      facet_to_coh_element->resize(new_nb_cohesive_elements);
      conn_cohesive->resize(new_nb_cohesive_elements);

      Array<Element> & event_list = element_event.getList();
      UInt event_list_old_size = event_list.getSize();
      event_list.resize(event_list_old_size + nb_facet_to_double);

      Element c_element(type_cohesive, 0, gt_facet, _ek_cohesive);
      Element f_element(type_facet, 0, gt_facet);

      UInt facets[2];

      for (UInt facet = 0; facet < nb_facet_to_double; ++facet) {

	/// (in 3D cohesive elements connectivity is inverted)
	facets[third_dimension] = f_to_double(facet);
	facets[!third_dimension] = old_nb_facet + facet;

	UInt cohesive_element = old_nb_cohesive_elements + facet;

	/// store doubled facets
	f_element.element = facets[0];
	(*facet_to_coh_element)(cohesive_element, 0) = f_element;
	f_element.element = facets[1];
	(*facet_to_coh_element)(cohesive_element, 1) = f_element;

	/// modify cohesive elements' connectivity
	for (UInt n = 0; n < nb_nodes_per_facet; ++n) {
	  (*conn_cohesive)(cohesive_element, n) = conn_facet(facets[0], n);
	  (*conn_cohesive)(cohesive_element, n + nb_nodes_per_facet)
	    = conn_facet(facets[1], n);
	}

	/// update element_to_facet vectors
	c_element.element = cohesive_element;
	element_to_facet(facets[0])[1] = c_element;
	element_to_facet(facets[1])[1] = c_element;

	/// add cohesive element to the element event list
	event_list(event_list_old_size + facet) = c_element;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::doublePointFacet(Mesh & mesh,
				 Mesh & mesh_facets,
				 Array<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  if (spatial_dimension != 1) return;

  Array<Real> & position = mesh.getNodes();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_facet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
    Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

    for(; it != end; ++it) {

      ElementType type_facet = *it;

      Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet, gt_facet);
      Array< std::vector<Element> > & element_to_facet
	= mesh_facets.getElementToSubelement(type_facet, gt_facet);

      const Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
								  type_facet,
								  gt_facet);

      UInt nb_facet_to_double = f_to_double.getSize();

      UInt old_nb_facet = element_to_facet.getSize() - nb_facet_to_double;
      UInt new_nb_facet = element_to_facet.getSize();

      UInt old_nb_nodes = position.getSize();
      UInt new_nb_nodes = old_nb_nodes + nb_facet_to_double;

      position.resize(new_nb_nodes);
      conn_facet.resize(new_nb_facet);

      UInt old_nb_doubled_nodes = doubled_nodes.getSize();
      doubled_nodes.resize(old_nb_doubled_nodes + nb_facet_to_double);

      for (UInt facet = 0; facet < nb_facet_to_double; ++facet) {
	UInt old_facet = f_to_double(facet);
	UInt new_facet = old_nb_facet + facet;

	ElementType type = element_to_facet(new_facet)[0].type;
	UInt el = element_to_facet(new_facet)[0].element;
	GhostType gt = element_to_facet(new_facet)[0].ghost_type;

	UInt old_node = conn_facet(old_facet);
	UInt new_node = old_nb_nodes + facet;

	/// update position
	position(new_node) = position(old_node);

	conn_facet(new_facet) = new_node;
	Array<UInt> & conn_segment = mesh.getConnectivity(type, gt);
	UInt nb_nodes_per_segment = conn_segment.getNbComponent();

	/// update facet connectivity
	UInt i;
	for (i = 0; conn_segment(el, i) != old_node
	       && i <= nb_nodes_per_segment; ++i);
	conn_segment(el, i) = new_node;

	doubled_nodes(old_nb_doubled_nodes + facet, 0) = old_node;
	doubled_nodes(old_nb_doubled_nodes + facet, 1) = new_node;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool third_dim_segments>
void MeshUtils::updateQuadraticSegments(Mesh & mesh,
					Mesh & mesh_facets,
					ElementType type_facet,
					GhostType gt_facet,
					Array<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  if (type_facet != _segment_3) return;

  Array<UInt> & f_to_double = mesh_facets.getData<UInt>("facet_to_double",
							type_facet,
							gt_facet);
  UInt nb_facet_to_double = f_to_double.getSize();

  UInt old_nb_facet
    = mesh_facets.getNbElement(type_facet, gt_facet) - nb_facet_to_double;

  Array<UInt> & conn_facet = mesh_facets.getConnectivity(type_facet,
							 gt_facet);

  Array< std::vector<Element> > & element_to_facet
    = mesh_facets.getElementToSubelement(type_facet, gt_facet);

  /// this ones matter only for segments in 3D
  Array< std::vector<Element> > * el_to_subfacet_double = NULL;
  Array< std::vector<Element> > * f_to_subfacet_double = NULL;

  if (third_dim_segments) {
    el_to_subfacet_double
      = & mesh_facets.getData<std::vector<Element> >("elements_to_subfacet_double",
						     type_facet,
						     gt_facet);

    f_to_subfacet_double
      = & mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
						     type_facet,
						     gt_facet);
  }

  Array<UInt> middle_nodes;
  middle_nodes.resize(0);

  for (UInt facet = 0; facet < nb_facet_to_double; ++facet) {
    UInt old_facet = f_to_double(facet);
    middle_nodes.push_back(conn_facet(old_facet, 2));
  }

  UInt old_nb_doubled_nodes = doubled_nodes.getSize();

  doubleNodes(mesh, middle_nodes, doubled_nodes);

  for (UInt n = old_nb_doubled_nodes; n < doubled_nodes.getSize(); ++n) {
    UInt old_node = doubled_nodes(n, 0);
    UInt new_node = doubled_nodes(n, 1);
    UInt sf = n - old_nb_doubled_nodes;
    UInt new_facet = old_nb_facet + sf;

    conn_facet(new_facet, 2) = new_node;

    if (third_dim_segments) {
      updateElementalConnectivity(mesh_facets,
				  old_node, new_node,
				  element_to_facet(new_facet));

      updateElementalConnectivity(mesh,
				  old_node, new_node,
				  (*el_to_subfacet_double)(sf),
				  &(*f_to_subfacet_double)(sf));
    }
    else {
      updateElementalConnectivity(mesh,
				  old_node, new_node,
				  element_to_facet(new_facet));
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::updateSubfacetToFacet(Mesh & mesh_facets,
				      ElementType type_subfacet,
				      GhostType gt_subfacet,
				      bool facet_mode) {
  AKANTU_DEBUG_IN();

  Array<UInt> & sf_to_double = mesh_facets.getData<UInt>("facet_to_double",
							 type_subfacet,
							 gt_subfacet);
  UInt nb_subfacet_to_double = sf_to_double.getSize();

  /// update subfacet_to_facet vector
  ElementType type_facet = _not_defined;
  GhostType gt_facet = _casper;
  Array<Element> * subfacet_to_facet = NULL;
  UInt nb_subfacet_per_facet = 0;
  UInt old_nb_subfacet
    = mesh_facets.getNbElement(type_subfacet, gt_subfacet) - nb_subfacet_to_double;

  Array<std::vector<Element> > * facet_list = NULL;
  if (facet_mode)
    facet_list = & mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
							      type_subfacet,
							      gt_subfacet);
  else
    facet_list = & mesh_facets.getData<std::vector<Element> >("subfacets_to_subsubfacet_double",
							      type_subfacet,
							      gt_subfacet);

  Element old_subfacet_el(type_subfacet, 0, gt_subfacet);
  Element new_subfacet_el(type_subfacet, 0, gt_subfacet);

  for (UInt sf = 0; sf < nb_subfacet_to_double; ++sf) {
    old_subfacet_el.element = sf_to_double(sf);
    new_subfacet_el.element = old_nb_subfacet + sf;

    for (UInt f = 0; f < (*facet_list)(sf).size(); ++f) {
      Element & facet = (*facet_list)(sf)[f];

      if (facet.type != type_facet || facet.ghost_type != gt_facet) {
	type_facet = facet.type;
	gt_facet = facet.ghost_type;

	subfacet_to_facet = & mesh_facets.getSubelementToElement(type_facet,
								 gt_facet);
	nb_subfacet_per_facet = subfacet_to_facet->getNbComponent();
      }

      Element * sf_update
	= std::find(subfacet_to_facet->storage()
		    + facet.element * nb_subfacet_per_facet,
		    subfacet_to_facet->storage()
		    + facet.element * nb_subfacet_per_facet + nb_subfacet_per_facet,
		    old_subfacet_el);

      AKANTU_DEBUG_ASSERT(subfacet_to_facet->storage()
			  + facet.element * nb_subfacet_per_facet !=
			  subfacet_to_facet->storage()
			  + facet.element * nb_subfacet_per_facet + nb_subfacet_per_facet,
			  "Subfacet not found");

      *sf_update = new_subfacet_el;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::updateFacetToSubfacet(Mesh & mesh_facets,
				      ElementType type_subfacet,
				      GhostType gt_subfacet,
				      bool facet_mode) {
  AKANTU_DEBUG_IN();

  Array<UInt> & sf_to_double = mesh_facets.getData<UInt>("facet_to_double",
							 type_subfacet,
							 gt_subfacet);
  UInt nb_subfacet_to_double = sf_to_double.getSize();

  Array< std::vector<Element> > & facet_to_subfacet
    = mesh_facets.getElementToSubelement(type_subfacet, gt_subfacet);

  Array< std::vector<Element> > * facet_to_subfacet_double = NULL;

  if (facet_mode) {
    facet_to_subfacet_double
      = & mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
						     type_subfacet,
						     gt_subfacet);
  }
  else {
    facet_to_subfacet_double
      = & mesh_facets.getData<std::vector<Element> >("subfacets_to_subsubfacet_double",
						     type_subfacet,
						     gt_subfacet);
  }

  for (UInt sf = 0; sf < nb_subfacet_to_double; ++sf)
    facet_to_subfacet.push_back((*facet_to_subfacet_double)(sf));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MeshUtils::doubleSubfacet(Mesh & mesh,
			       Mesh & mesh_facets,
			       Array<UInt> & doubled_nodes) {
  AKANTU_DEBUG_IN();

  if (spatial_dimension == 1) return;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType gt_subfacet = *gt;

    Mesh::type_iterator it  = mesh_facets.firstType(0, gt_subfacet);
    Mesh::type_iterator end = mesh_facets.lastType(0, gt_subfacet);

    for(; it != end; ++it) {

      ElementType type_subfacet = *it;

      Array<UInt> & sf_to_double = mesh_facets.getData<UInt>("facet_to_double",
							     type_subfacet,
							     gt_subfacet);
      UInt nb_subfacet_to_double = sf_to_double.getSize();

      if (nb_subfacet_to_double == 0) continue;

      AKANTU_DEBUG_ASSERT(type_subfacet == _point_1,
			  "Only _point_1 subfacet doubling is supported at the moment");

      Array<UInt> & conn_subfacet = mesh_facets.getConnectivity(type_subfacet,
								gt_subfacet);

      UInt old_nb_subfacet = conn_subfacet.getSize();
      UInt new_nb_subfacet = old_nb_subfacet + nb_subfacet_to_double;

      conn_subfacet.resize(new_nb_subfacet);

      Array<UInt> nodes_to_double;
      UInt old_nb_doubled_nodes = doubled_nodes.getSize();

      /// double nodes
      for (UInt sf = 0; sf < nb_subfacet_to_double; ++sf) {
	UInt old_subfacet = sf_to_double(sf);
	nodes_to_double.push_back(conn_subfacet(old_subfacet));
      }

      doubleNodes(mesh, nodes_to_double, doubled_nodes);

      /// add new nodes in connectivity
      for (UInt sf = 0; sf < nb_subfacet_to_double; ++sf) {
	UInt new_subfacet = old_nb_subfacet + sf;
	UInt new_node = doubled_nodes(old_nb_doubled_nodes + sf, 1);

	conn_subfacet(new_subfacet) = new_node;
      }

      /// update facet and element connectivity
      Array<std::vector<Element> > & f_to_subfacet_double =
	mesh_facets.getData<std::vector<Element> >("facets_to_subfacet_double",
						   type_subfacet,
						   gt_subfacet);

      Array<std::vector<Element> > & el_to_subfacet_double =
	mesh_facets.getData<std::vector<Element> >("elements_to_subfacet_double",
						   type_subfacet,
						   gt_subfacet);

      Array<std::vector<Element> > * sf_to_subfacet_double = NULL;

      if (spatial_dimension == 3)
	sf_to_subfacet_double
	  = & mesh_facets.getData<std::vector<Element> >("subfacets_to_subsubfacet_double",
							 type_subfacet,
							 gt_subfacet);

      for (UInt sf = 0; sf < nb_subfacet_to_double; ++sf) {
	UInt old_node = doubled_nodes(old_nb_doubled_nodes + sf, 0);
	UInt new_node = doubled_nodes(old_nb_doubled_nodes + sf, 1);

	updateElementalConnectivity(mesh,
				    old_node, new_node,
				    el_to_subfacet_double(sf),
				    &f_to_subfacet_double(sf));

	updateElementalConnectivity(mesh_facets,
				    old_node, new_node,
				    f_to_subfacet_double(sf));

	if (spatial_dimension == 3)
	  updateElementalConnectivity(mesh_facets,
				      old_node, new_node,
				      (*sf_to_subfacet_double)(sf));
      }

      if (spatial_dimension == 2) {
	updateSubfacetToFacet(mesh_facets, type_subfacet, gt_subfacet, true);
	updateFacetToSubfacet(mesh_facets, type_subfacet, gt_subfacet, true);
      }
      else if (spatial_dimension == 3) {
	updateSubfacetToFacet(mesh_facets, type_subfacet, gt_subfacet, false);
	updateFacetToSubfacet(mesh_facets, type_subfacet, gt_subfacet, false);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::flipFacets(Mesh & mesh_facets,
			   const ByElementTypeUInt & global_connectivity,
			   GhostType gt_facet) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();

  Mesh::type_iterator it  = mesh_facets.firstType(spatial_dimension - 1, gt_facet);
  Mesh::type_iterator end = mesh_facets.lastType(spatial_dimension - 1, gt_facet);

  /// loop on every ghost facet
  for(; it != end; ++it) {
    ElementType type_facet = *it;

    Array<UInt> & connectivity = mesh_facets.getConnectivity(type_facet, gt_facet);
    const Array<UInt> & g_connectivity = global_connectivity(type_facet, gt_facet);

    Array<std::vector<Element> > & el_to_f =
      mesh_facets.getElementToSubelement(type_facet, gt_facet);
    Array<Element> & subfacet_to_facet =
      mesh_facets.getSubelementToElement(type_facet, gt_facet);

    UInt nb_subfacet_per_facet = subfacet_to_facet.getNbComponent();
    UInt nb_nodes_per_facet = connectivity.getNbComponent();
    UInt nb_facet = connectivity.getSize();

    UInt nb_nodes_per_P1_facet
      = Mesh::getNbNodesPerElement(Mesh::getP1ElementType(type_facet));

    /// get global connectivity for local mesh
    Array<UInt> global_conn_tmp(nb_facet, nb_nodes_per_facet);
    mesh_facets.getGlobalConnectivity(global_conn_tmp, type_facet, gt_facet);

    Array<UInt>::iterator<Vector<UInt> > conn_it =
      connectivity.begin(nb_nodes_per_facet);
    Array<UInt>::iterator<Vector<UInt> > gconn_tmp_it =
      global_conn_tmp.begin(nb_nodes_per_facet);
    Array<UInt>::const_iterator<Vector<UInt> > conn_glob_it =
      g_connectivity.begin(nb_nodes_per_facet);
    Array<Element>::iterator<Vector<Element> > subf_to_f =
      subfacet_to_facet.begin(nb_subfacet_per_facet);

    UInt * conn_tmp_pointer = new UInt[nb_nodes_per_facet];
    Vector<UInt> conn_tmp(conn_tmp_pointer, nb_nodes_per_facet);

    Element el_tmp;
    Element * subf_tmp_pointer = new Element[nb_subfacet_per_facet];
    Vector<Element> subf_tmp(subf_tmp_pointer, nb_subfacet_per_facet);

    for (UInt f = 0; f < nb_facet; ++f, ++conn_it, ++gconn_tmp_it,
	   ++subf_to_f, ++conn_glob_it) {

      Vector<UInt> & gconn_tmp  = *gconn_tmp_it;
      const Vector<UInt> & conn_glob  = *conn_glob_it;
      Vector<UInt> & conn_local = *conn_it;

      /// skip facet if connectivities are the same
      if (gconn_tmp == conn_glob) continue;

      /// re-arrange connectivity
      conn_tmp = conn_local;

      UInt * begin = conn_glob.storage();
      UInt * end = conn_glob.storage() + nb_nodes_per_facet;

      for (UInt n = 0; n < nb_nodes_per_facet; ++n) {

	UInt * new_node = std::find(begin, end, gconn_tmp(n));
	AKANTU_DEBUG_ASSERT(new_node != end, "Node not found");

	UInt new_position = new_node - begin;

	conn_local(new_position) = conn_tmp(n);
      }

      /// if 3D, check if facets are just rotated
      if (spatial_dimension == 3) {
	/// find first node
	UInt * new_node = std::find(begin, end, gconn_tmp(0));
	AKANTU_DEBUG_ASSERT(new_node != end, "Node not found");

	UInt new_position = new_node - begin;
	UInt n = 1;

	/// count how many nodes in the received connectivity follow
	/// the same order of those in the local connectivity
	for (; n < nb_nodes_per_P1_facet &&
	       gconn_tmp(n) == conn_glob((new_position + n) % nb_nodes_per_P1_facet);
	     ++n);

	/// skip the facet inversion if facet is just rotated
	if (n == nb_nodes_per_P1_facet) continue;
      }

      /// update data to invert facet
      el_tmp = el_to_f(f)[0];
      el_to_f(f)[0] = el_to_f(f)[1];
      el_to_f(f)[1] = el_tmp;

      subf_tmp = (*subf_to_f);
      (*subf_to_f)(0) = subf_tmp(1);
      (*subf_to_f)(1) = subf_tmp(0);
    }

    delete [] subf_tmp_pointer;
    delete [] conn_tmp_pointer;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MeshUtils::fillElementToSubElementsData(Mesh & mesh) {
  AKANTU_DEBUG_IN();

  if(mesh.getNbElement(mesh.getSpatialDimension() - 1) == 0) {
    AKANTU_DEBUG_WARNING("There are not facets, add them in the mesh file or call the buildFacet method.");
    return;
  }

  UInt spatial_dimension = mesh.getSpatialDimension();
  ByElementTypeReal barycenters("barycenter_tmp", mesh.getID());
  mesh.initByElementTypeArray(barycenters,
                              spatial_dimension,
                              _all_dimensions);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end();
       ++gt) {
    Mesh::type_iterator it  = mesh.firstType(_all_dimensions, *gt);
    Mesh::type_iterator end = mesh.lastType(_all_dimensions, *gt);
    for(; it != end; ++it) {
      UInt nb_element = mesh.getNbElement(*it, *gt);
      Array<Real> & barycenters_arr = barycenters(*it, *gt);
      barycenters_arr.resize(nb_element);
      Array<Real>::iterator< Vector<Real> > bary = barycenters_arr.begin(spatial_dimension);
      Array<Real>::iterator< Vector<Real> > bary_end = barycenters_arr.end(spatial_dimension);

      for (UInt el = 0; bary != bary_end; ++bary, ++el) {
        mesh.getBarycenter(el, *it, bary->storage(), *gt);
      }
    }
  }

  for(Int sp(spatial_dimension); sp >= 1; --sp) {
    if(mesh.getNbElement(sp) == 0) continue;

    for (ghost_type_t::iterator git = ghost_type_t::begin();  git != ghost_type_t::end(); ++git) {
      Mesh::type_iterator tit  = mesh.firstType(sp, *git);
      Mesh::type_iterator tend = mesh.lastType(sp, *git);
      for (;tit != tend; ++tit) {
        mesh.getSubelementToElementPointer(*tit, *git)->resize(mesh.getNbElement(*tit, *git));
        mesh.getSubelementToElement(*tit, *git).clear();
      }

      tit  = mesh.firstType(sp - 1, *git);
      tend = mesh.lastType(sp - 1, *git);
      for (;tit != tend; ++tit) {
        mesh.getElementToSubelementPointer(*tit, *git)->resize(mesh.getNbElement(*tit, *git));
        mesh.getElementToSubelement(*tit, *git).clear();
      }
    }


    CSR<Element> nodes_to_elements;
    buildNode2Elements(mesh, nodes_to_elements, sp);

    Element facet_element;

    for (ghost_type_t::iterator git = ghost_type_t::begin();  git != ghost_type_t::end(); ++git) {
      Mesh::type_iterator tit  = mesh.firstType(sp - 1, *git);
      Mesh::type_iterator tend = mesh.lastType(sp - 1, *git);

      facet_element.ghost_type = *git;
      for (;tit != tend; ++tit) {
        facet_element.type = *tit;

        Array< std::vector<Element> > & element_to_subelement = mesh.getElementToSubelement(*tit, *git);

        const Array<UInt> & connectivity = mesh.getConnectivity(*tit, *git);

        Array<UInt>::const_iterator< Vector<UInt> > fit  = connectivity.begin(mesh.getNbNodesPerElement(*tit));
        Array<UInt>::const_iterator< Vector<UInt> > fend = connectivity.end(mesh.getNbNodesPerElement(*tit));

        UInt fid = 0;
        for (;fit != fend; ++fit, ++fid) {
          const Vector<UInt> & facet = *fit;
          facet_element.element = fid;
          std::map<Element, UInt> element_seen_counter;
          UInt nb_nodes_per_facet = mesh.getNbNodesPerElement(Mesh::getP1ElementType(*tit));
          for (UInt n(0); n < nb_nodes_per_facet; ++n) {
            CSR<Element>::iterator eit  = nodes_to_elements.begin(facet(n));
            CSR<Element>::iterator eend = nodes_to_elements.end(facet(n));
            for(;eit != eend; ++eit) {
              Element & elem = *eit;
              std::map<Element, UInt>::iterator cit = element_seen_counter.find(elem);
              if(cit != element_seen_counter.end()) {
                cit->second++;
              } else {
                element_seen_counter[elem] = 1;
              }
            }
          }

          std::vector<Element> connected_elements;
          std::map<Element, UInt>::iterator cit  = element_seen_counter.begin();
          std::map<Element, UInt>::iterator cend = element_seen_counter.end();
          for(;cit != cend; ++cit) {
            if(cit->second == nb_nodes_per_facet) connected_elements.push_back(cit->first);
          }
          if(connected_elements.size() == 1)
            MeshUtils::sortElements(connected_elements, facet, mesh, mesh, barycenters);

          std::vector<Element>::iterator ceit  = connected_elements.begin();
          std::vector<Element>::iterator ceend = connected_elements.end();
          for(;ceit != ceend; ++ceit)
            element_to_subelement(fid).push_back(*ceit);

          for (UInt ce = 0; ce < connected_elements.size(); ++ce) {
            Element & elem = connected_elements[ce];
            Array<Element> & subelement_to_element =
              *(mesh.getSubelementToElementPointer(elem.type,
                                                   elem.ghost_type));

            UInt f(0);
            for(; f < mesh.getNbFacetsPerElement(elem.type) && subelement_to_element(elem.element, f) != ElementNull; ++f);
	    AKANTU_DEBUG_ASSERT(f < mesh.getNbFacetsPerElement(elem.type), "The element "<< elem << " seems to have too many facets!!");
            subelement_to_element(elem.element, f) = facet_element;
          }
        }
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <bool third_dim_points>
bool MeshUtils::findElementsAroundSubfacet(const Mesh & mesh,
					   const Mesh & mesh_facets,
					   const Element & starting_element,
					   const Element & end_facet,
					   const Vector<UInt> & subfacet_connectivity,
					   std::vector<Element> & elem_list,
					   std::vector<Element> & facet_list,
					   std::vector<Element> * subfacet_list) {
  AKANTU_DEBUG_IN();

  /// preallocated stuff before starting
  bool facet_matched = false;

  elem_list.resize(0);
  facet_list.resize(0);

  if (third_dim_points)
    subfacet_list->resize(0);

  elem_list.push_back(starting_element);

  const Array<UInt> * facet_connectivity = NULL;
  const Array<UInt> * sf_connectivity = NULL;
  const Array<Element> * facet_to_element = NULL;
  const Array<Element> * subfacet_to_facet = NULL;

  ElementType current_type = _not_defined;
  GhostType current_ghost_type = _casper;

  ElementType current_facet_type = _not_defined;
  GhostType current_facet_ghost_type = _casper;

  ElementType current_subfacet_type = _not_defined;
  GhostType current_subfacet_ghost_type = _casper;

  const Array< std::vector<Element> > * element_to_facet = NULL;

  const Element * opposing_el = NULL;

  Array<Element> elements_to_check;
  elements_to_check.resize(0);
  elements_to_check.push_back(starting_element);

  /// keep going until there are elements to check
  while (elements_to_check.getSize() != 0) {

    /// loop over each element to check
    for (UInt el = 0; el < elements_to_check.getSize(); ++el) {

      Element & current_el = elements_to_check(el);

      if (current_el.type != current_type ||
	  current_el.ghost_type != current_ghost_type) {

	current_type = current_el.type;
	current_ghost_type = current_el.ghost_type;

	facet_to_element = & mesh_facets.getSubelementToElement(current_type,
								current_ghost_type);
      }

      /// loop over each facet of the element
      for (UInt f = 0; f < facet_to_element->getNbComponent(); ++f) {

	const Element & current_facet = (*facet_to_element)(current_el.element, f);

	if (current_facet == ElementNull) continue;

	if (current_facet_type != current_facet.type ||
	    current_facet_ghost_type != current_facet.ghost_type) {

	  current_facet_type = current_facet.type;
	  current_facet_ghost_type = current_facet.ghost_type;

	  element_to_facet =
	    & mesh_facets.getElementToSubelement(current_facet_type,
						 current_facet_ghost_type);
	  facet_connectivity = & mesh_facets.getConnectivity(current_facet_type,
							     current_facet_ghost_type);

	  if (third_dim_points)
	    subfacet_to_facet =
	      & mesh_facets.getSubelementToElement(current_facet_type,
						   current_facet_ghost_type);
	}

	/// check if end facet is reached
	if (current_facet == end_facet)
	  facet_matched = true;

	/// add this facet if not already passed
	if (std::find(facet_list.begin(),
		      facet_list.end(),
		      current_facet) == facet_list.end() &&
	    hasElement(*facet_connectivity, current_facet, subfacet_connectivity)) {
	  facet_list.push_back(current_facet);

	  if (third_dim_points) {
	    /// check subfacets
	    for (UInt sf = 0; sf < subfacet_to_facet->getNbComponent(); ++sf) {
	      const Element & current_subfacet
		= (*subfacet_to_facet)(current_facet.element, sf);

	      if (current_subfacet == ElementNull) continue;

	      if (current_subfacet_type != current_subfacet.type ||
		  current_subfacet_ghost_type != current_subfacet.ghost_type) {
		current_subfacet_type = current_subfacet.type;
		current_subfacet_ghost_type = current_subfacet.ghost_type;

		sf_connectivity
		  = & mesh_facets.getConnectivity(current_subfacet_type,
						  current_subfacet_ghost_type);
	      }

	      if (std::find(subfacet_list->begin(),
			    subfacet_list->end(),
			    current_subfacet) == subfacet_list->end() &&
		  hasElement(*sf_connectivity, current_subfacet, subfacet_connectivity))
		subfacet_list->push_back(current_subfacet);
	    }
	  }
	}
	else
	  continue;

	/// consider opposing element
	if ( (*element_to_facet)(current_facet.element)[0] == current_el)
	  opposing_el = & (*element_to_facet)(current_facet.element)[1];
	else
	  opposing_el = & (*element_to_facet)(current_facet.element)[0];

	/// skip null elements since they are on a boundary
	if (*opposing_el == ElementNull) continue;

	/// skip this element if already added
	if ( std::find(elem_list.begin(),
		       elem_list.end(),
		       *opposing_el) != elem_list.end() ) continue;

	/// only regular elements have to be checked
	if (opposing_el->kind == _ek_regular)
	  elements_to_check.push_back(*opposing_el);

	elem_list.push_back(*opposing_el);

#ifndef AKANTU_NDEBUG
	const Array<UInt> & conn_elem = mesh.getConnectivity(opposing_el->type,
							     opposing_el->ghost_type);

	AKANTU_DEBUG_ASSERT(hasElement(conn_elem, *opposing_el, subfacet_connectivity),
			    "Subfacet doesn't belong to this element");
#endif
      }

      /// erased checked element from the list
      elements_to_check.erase(el);
    }
  }

  AKANTU_DEBUG_OUT();
  return facet_matched;
}

/* -------------------------------------------------------------------------- */
/* Internal functions                                                         */
/* -------------------------------------------------------------------------- */

void MeshUtils::sortElements(std::vector<Element> & elements, const Vector<UInt> facet,
                             const Mesh & mesh, const Mesh & mesh_facets,
                             const ByElementTypeReal & barycenters) {
  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_element_connected_to_facet = elements.size();

  const Array<Real> & mesh_facets_nodes = mesh_facets.getNodes();
  const Array<Real>::const_iterator< Vector<Real> > mesh_facets_nodes_it =
    mesh_facets_nodes.begin(spatial_dimension);

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
      barycenters(elements[0].type, elements[0].ghost_type).begin(spatial_dimension);

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
        barycenters(elements[n].type, elements[n].ghost_type).begin(spatial_dimension);
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
        barycenters(elements[n].type, elements[n].ghost_type).begin(spatial_dimension);
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

__END_AKANTU__

//  LocalWords:  ElementType
