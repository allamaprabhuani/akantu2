/**
 * @file   contact_search_3d_explicit.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:49:04 2010
 *
 * @brief  Specialization of  the  contact  search structure  for  3D within  an
 * explicit time scheme
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
#include "regular_grid_neighbor_structure.hh"
#include "contact_search_3d_explicit.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactSearch3dExplicit::ContactSearch3dExplicit(Contact & contact,
						 const ContactNeighborStructureType & neighbors_structure_type,
						 const ContactSearchType & type,
						 const ContactSearchID & id) :
  ContactSearch(contact, neighbors_structure_type, type, id), spatial_dimension(contact.getModel().getSpatialDimension()), mesh(contact.getModel().getFEM().getMesh()) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::findPenetration(const Surface & master_surface, PenetrationList & penetration_list) {
  AKANTU_DEBUG_IN();

  /// get the NodesNeighborList for the given master surface
  std::map<Surface, ContactNeighborStructure *>::iterator it_surface;
  for (it_surface = neighbors_structure.begin(); it_surface != neighbors_structure.end(); ++it_surface) {
    if(it_surface->first == master_surface) {
      break;
    }
  }
  AKANTU_DEBUG_ASSERT(it_surface != neighbors_structure.end(),
		      "Master surface not found in this search object " << id);
  const NodesNeighborList & neighbor_list = dynamic_cast<const NodesNeighborList&>(it_surface->second->getNeighborList());
  UInt nb_impactor_nodes = neighbor_list.impactor_nodes.getSize();

  Vector<UInt> * closest_master_nodes = new Vector<UInt>(nb_impactor_nodes, 1);
  Vector<bool> * has_closest_master_node = new Vector<bool>(nb_impactor_nodes, 1, false);
  findClosestMasterNodes(master_surface, closest_master_nodes, has_closest_master_node);
  UInt * closest_master_nodes_val = closest_master_nodes->values;
  bool * has_closest_master_node_val = has_closest_master_node->values;

  /// get list of impactor and master nodes from neighbor list
  UInt * impactor_nodes_val = neighbor_list.impactor_nodes.values;

  //const Mesh & mesh = contact.getModel().getFEM().getMesh();
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  /// find existing surface element types
  UInt nb_types = type_list.size();
  UInt nb_facet_types = 0;
  ElementType facet_type[_max_element_type];

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

    if(mesh.getSpatialDimension(type) == spatial_dimension) {

      ElementType current_facet_type = mesh.getFacetElementType(type);
      facet_type[nb_facet_types++] = current_facet_type;

      /// initialization of penetration list
      std::stringstream sstr_facets_offset;
      sstr_facets_offset << id << ":penetrated_facets_offset:" << current_facet_type;
      penetration_list.penetrated_facets_offset[current_facet_type] = new Vector<UInt>(0, 1, sstr_facets_offset.str());

      std::stringstream sstr_facets;
      sstr_facets << id << ":penetrated_facets:" << current_facet_type;
      penetration_list.penetrated_facets[current_facet_type] = new Vector<UInt>(0, 1, sstr_facets.str());

      std::stringstream sstr_normals;
      sstr_normals << id << ":facets_normals:" << current_facet_type;
      penetration_list.facets_normals[current_facet_type] = new Vector<Real>(0, spatial_dimension, sstr_normals.str());

      std::stringstream sstr_gaps;
      sstr_gaps << id << ":gaps:" << current_facet_type;
      penetration_list.gaps[current_facet_type] = new Vector<Real>(0, 1, sstr_gaps.str());

      std::stringstream sstr_projected_positions;
      sstr_projected_positions << id << ":projected_positions:" << current_facet_type;
      penetration_list.projected_positions[current_facet_type] = new Vector<Real>(0, spatial_dimension, sstr_projected_positions.str());
    }
  }

  for(UInt in = 0; in < nb_impactor_nodes; ++in) {

    if (!has_closest_master_node_val[in])
      continue;

    UInt current_impactor_node = impactor_nodes_val[in];
    UInt closest_master_node = closest_master_nodes_val[in];

    std::vector<Element> surface_elements;
    Vector<bool> * are_inside = new Vector<bool>(0, 1);
    Vector<bool> * are_in_projection_area = new Vector<bool>(0, 1);

    Element considered_element;

    for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {
      ElementType type = facet_type[el_type];

      UInt * surface_id_val = mesh.getSurfaceId(type).values;

      const Vector<UInt> & node_to_elements_offset = contact.getNodeToElementsOffset(type);
      const Vector<UInt> & node_to_elements = contact.getNodeToElements(type);
      UInt * node_to_elements_offset_val = node_to_elements_offset.values;
      UInt * node_to_elements_val        = node_to_elements.values;

      UInt min_element_offset = node_to_elements_offset_val[closest_master_node];
      UInt max_element_offset = node_to_elements_offset_val[closest_master_node + 1];

      considered_element.type = type;

      for(UInt el = min_element_offset; el < max_element_offset; ++el) {
	UInt surface_element = node_to_elements_val[el];
	if(surface_id_val[surface_element] == master_surface) {
	  bool is_inside;
	  bool is_in_projection_area;
	  checkPenetrationSituation(current_impactor_node,
				    surface_element,
				    type,
				    is_inside,
				    is_in_projection_area);

	  considered_element.element = surface_element;

	  surface_elements.push_back(considered_element);
	  are_inside->push_back(is_inside);
	  are_in_projection_area->push_back(is_in_projection_area);
	}
      }
    }

    UInt nb_penetrated_elements = 0;
    UInt nb_elements_type[_max_element_type]; memset(nb_elements_type, 0, sizeof(UInt) * _max_element_type);

    bool * are_inside_val = are_inside->values;
    bool * are_in_projection_area_val = are_in_projection_area->values;

    bool is_outside = false;
    UInt nb_surface_elements = static_cast<UInt>(surface_elements.size());

    // add all elements for which impactor is inside and in projection area
    for(UInt el = 0; el < nb_surface_elements; ++el) {
      if(are_inside_val[el] && are_in_projection_area_val[el]) {

	ElementType current_type = surface_elements.at(el).type;
	UInt current_element = surface_elements.at(el).element;
	penetration_list.penetrated_facets[current_type]->push_back(current_element);

	Real normal[3];
	Real projected_position[3];
	Real gap;
	computeComponentsOfProjection(current_impactor_node,
				      current_element,
				      current_type,
				      normal,
				      gap,
				      projected_position);

	penetration_list.facets_normals[current_type]->push_back(normal);
	penetration_list.projected_positions[current_type]->push_back(projected_position);
	penetration_list.gaps[current_type]->push_back(gap);

	nb_penetrated_elements++;
	nb_elements_type[current_type]++;
      }
    }
    if(nb_penetrated_elements > 0) {
      for(it = type_list.begin(); it != type_list.end(); ++it) {
	ElementType type = *it;
	if(mesh.getSpatialDimension(type) == spatial_dimension) {
	  penetration_list.penetrating_nodes.push_back(current_impactor_node);
	  ElementType current_facet_type = mesh.getFacetElementType(type);
	  penetration_list.penetrated_facets_offset[current_facet_type]->push_back(nb_elements_type[current_facet_type]);
	}
      }
    }

    // if there was no element which is inside and in projection area
    // check if node is not definitly outside
    else {
      for(UInt el = 0; el < nb_surface_elements && !is_outside; ++el) {
	if(!are_inside_val[el] && are_in_projection_area_val[el]) {
	  is_outside = true;
	}
      }

      // it is not definitly outside take all elements to which it is at least inside
      if(!is_outside) {
	bool found_inside_node = false;
	for(UInt el = 0; el < nb_surface_elements; ++el) {
	  if(are_inside_val[el] && !are_in_projection_area_val[el]) {

	    ElementType current_type = surface_elements.at(el).type;
	    UInt current_element = surface_elements.at(el).element;
	    penetration_list.penetrated_facets[current_type]->push_back(current_element);

	    Real normal[3];
	    Real projected_position[3];
	    Real gap;
	    computeComponentsOfProjection(current_impactor_node,
					  current_element,
					  current_type,
					  normal,
					  gap,
					  projected_position);

	    penetration_list.facets_normals[current_type]->push_back(normal);
	    penetration_list.projected_positions[current_type]->push_back(projected_position);
	    penetration_list.gaps[current_type]->push_back(gap);

	    nb_penetrated_elements++;
	    nb_elements_type[current_type]++;
	    found_inside_node = true;
	  }
	}
	if(found_inside_node) {
	  for(it = type_list.begin(); it != type_list.end(); ++it) {
	    ElementType type = *it;
	    if(mesh.getSpatialDimension(type) == spatial_dimension) {
	      penetration_list.penetrating_nodes.push_back(current_impactor_node);
	      ElementType current_facet_type = mesh.getFacetElementType(type);
	      penetration_list.penetrated_facets_offset[current_facet_type]->push_back(nb_elements_type[current_facet_type]);
	    }
	  }
	}
      }
    }

    delete are_in_projection_area;
    delete are_inside;
  }

  // make the offset table a real offset table
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;

    if(mesh.getSpatialDimension(type) == spatial_dimension) {
      ElementType current_facet_type = mesh.getFacetElementType(type);

      UInt tmp_nb_facets = penetration_list.penetrated_facets_offset[current_facet_type]->getSize();
      penetration_list.penetrated_facets_offset[current_facet_type]->resize(tmp_nb_facets+1);

      Vector<UInt> & tmp_penetrated_facets_offset = *(penetration_list.penetrated_facets_offset[current_facet_type]);
      UInt * tmp_penetrated_facets_offset_val = tmp_penetrated_facets_offset.values;

      for (UInt i = 1; i < tmp_nb_facets; ++i)
	tmp_penetrated_facets_offset_val[i] += tmp_penetrated_facets_offset_val[i - 1];
      for (UInt i = tmp_nb_facets; i > 0; --i)
	tmp_penetrated_facets_offset_val[i] = tmp_penetrated_facets_offset_val[i - 1];
      tmp_penetrated_facets_offset_val[0] = 0;
    }
  }

  delete closest_master_nodes;
  delete has_closest_master_node;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::findClosestMasterNodes(const Surface & master_surface,
						     Vector<UInt> * closest_master_nodes,
						     Vector<bool> * has_closest_master_node) {
  AKANTU_DEBUG_IN();

  bool * has_closest_master_node_val = has_closest_master_node->values;
  UInt * closest_master_nodes_val = closest_master_nodes->values;

  /// get the NodesNeighborList for the given master surface
  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    if(it->first == master_surface) {
      break;
    }
  }
  AKANTU_DEBUG_ASSERT(it != neighbors_structure.end(), "Master surface not found in this search object " << id);
  const NodesNeighborList & neighbor_list = dynamic_cast<const NodesNeighborList&>(it->second->getNeighborList());

  /// get list of impactor and master nodes from neighbor list
  UInt * impactor_nodes_val = neighbor_list.impactor_nodes.values;
  UInt * master_nodes_offset_val = neighbor_list.master_nodes_offset.values;
  UInt * master_nodes_val = neighbor_list.master_nodes.values;

  /// loop over all impactor nodes and find for each the closest master node
  UInt nb_impactor_nodes = neighbor_list.impactor_nodes.getSize();
  for(UInt imp = 0; imp < nb_impactor_nodes; ++imp) {
    UInt current_impactor_node = impactor_nodes_val[imp];
    UInt min_offset = master_nodes_offset_val[imp];
    UInt max_offset = master_nodes_offset_val[imp + 1];

    // if there is no master node go to next impactor node
    if (min_offset == max_offset)
      continue;

    Real min_square_distance = std::numeric_limits<Real>::max();
    UInt closest_master_node = (UInt)-1;                     // for finding error
    for(UInt mn = min_offset; mn < max_offset; ++mn) {
      UInt current_master_node = master_nodes_val[mn];
      Real square_distance = computeSquareDistanceBetweenNodes(current_impactor_node, current_master_node);
      if(min_square_distance > square_distance) {
	min_square_distance = square_distance;
	closest_master_node = current_master_node;
      }
    }
    AKANTU_DEBUG_ASSERT(closest_master_node != ((UInt)-1), "could not find a closest master node for impactor node: " << current_impactor_node);
    has_closest_master_node_val[imp] = true;
    closest_master_nodes_val[imp] = closest_master_node;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::computeComponentsOfProjection(const UInt impactor_node,
							    const UInt surface_element,
							    const ElementType type,
							    Real * normal,
							    Real & gap,
							    Real * projected_position) {
  AKANTU_DEBUG_IN();

  switch(type) {
  case _segment_2: {
    computeComponentsOfProjectionSegment2(impactor_node, surface_element, normal, gap, projected_position);
    break;
  }
  case _triangle_3: {
    computeComponentsOfProjectionTriangle3(impactor_node, surface_element, normal, gap, projected_position);
    break;
  }
  case _not_defined: {
    AKANTU_DEBUG_ERROR("Not a valid surface element type : " << type);
    break;
  }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::checkPenetrationSituation(const UInt impactor_node,
							const UInt surface_element,
							const ElementType type,
							bool & is_inside,
							bool & is_in_projection_area) {
  AKANTU_DEBUG_IN();

  switch(type) {
  case _segment_2: {
    checkPenetrationSituationSegment2(impactor_node, surface_element, is_inside, is_in_projection_area);
    break;
  }
  case _triangle_3: {
    checkPenetrationSituationTriangle3(impactor_node, surface_element, is_inside, is_in_projection_area);
    break;
  }
  case _not_defined: {
    AKANTU_DEBUG_ERROR("Not a valid surface element type : " << type);
    break;
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::computeComponentsOfProjectionSegment2(const UInt impactor_node,
								    const UInt surface_element,
								    Real * normal,
								    Real & gap,
								    Real * projected_position) {
  AKANTU_DEBUG_IN();

  const UInt dim = spatial_dimension;
  const ElementType type = _segment_2;
  const UInt nb_nodes_element = Mesh::getNbNodesPerElement(type);

  Real * current_position = contact.getModel().getCurrentPosition().values;
  UInt * connectivity = mesh.getConnectivity(type).values;

  UInt node_1 = surface_element * nb_nodes_element;
  Real * position_node_1 = &(current_position[connectivity[node_1 + 0] * spatial_dimension]);
  Real * position_node_2 = &(current_position[connectivity[node_1 + 1] * spatial_dimension]);
  Real * position_impactor_node = &(current_position[impactor_node * spatial_dimension]);

  /// compute the normal of the face
  Real vector_1[2];
  Math::vector_2d(position_node_1, position_node_2, vector_1); /// @todo: check if correct order of nodes !!
  Math::normal2(vector_1, normal);

  /// compute the gap between impactor and face
  /// gap positive if impactor outside of surface
  Real vector_node_1_impactor[2];
  Math::vector_2d(position_node_1, position_impactor_node, vector_node_1_impactor);
  gap = Math::vectorDot2(vector_node_1_impactor, normal);

  /// compute the projected position of the impactor node onto the face
  for(UInt i=0; i < dim; ++i) {
    projected_position[i] = position_impactor_node[i] - gap * normal[i];
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::computeComponentsOfProjectionTriangle3(const UInt impactor_node,
								     const UInt surface_element,
								     Real * normal,
								     Real & gap,
								     Real * projected_position) {
  AKANTU_DEBUG_IN();

  const UInt dim = spatial_dimension;
  const ElementType type = _triangle_3;
  const UInt nb_nodes_element = Mesh::getNbNodesPerElement(type);

  Real * current_position = contact.getModel().getCurrentPosition().values;
  UInt * connectivity = mesh.getConnectivity(type).values;

  UInt node_1 = surface_element * nb_nodes_element;
  Real * position_node_1 = &(current_position[connectivity[node_1 + 0] * spatial_dimension]);
  Real * position_node_2 = &(current_position[connectivity[node_1 + 1] * spatial_dimension]);
  Real * position_node_3 = &(current_position[connectivity[node_1 + 2] * spatial_dimension]);
  Real * position_impactor_node = &(current_position[impactor_node * spatial_dimension]);

  /// compute the normal of the face
  Real vector_1[3];
  Real vector_2[3];
  Math::vector_3d(position_node_1, position_node_2, vector_1); /// @todo: check if correct order of nodes !!
  Math::vector_3d(position_node_1, position_node_3, vector_2);
  Math::normal3(vector_1, vector_2, normal);

  /// compute the gap between impactor and face
  /// gap positive if impactor outside of surface
  Real vector_node_1_impactor[3];
  Math::vector_3d(position_node_1, position_impactor_node, vector_node_1_impactor);
  gap = Math::vectorDot3(vector_node_1_impactor, normal);

  /// compute the projected position of the impactor node onto the face
  for(UInt i=0; i < dim; ++i) {
    projected_position[i] = position_impactor_node[i] - gap * normal[i];
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::checkPenetrationSituationSegment2(const UInt impactor_node,
								const UInt surface_element,
								bool & is_inside,
								bool & is_in_projection_area) {

  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(spatial_dimension == 2, "wrong spatial dimension (=" << spatial_dimension << ") for checkPenetrationSituationSegment2");
  const UInt dim = spatial_dimension;
  const ElementType type = _segment_2;
  const UInt nb_nodes_element = Mesh::getNbNodesPerElement(type);
  const Real tolerance = std::numeric_limits<Real>::epsilon();

  Real * current_position = contact.getModel().getCurrentPosition().values;
  UInt * connectivity = mesh.getConnectivity(type).values;

  Real gap;
  Real normal[2];
  Real projected_position[2];
  computeComponentsOfProjectionSegment2(impactor_node, surface_element, normal, gap, projected_position);

  // -------------------------------------------------------
  /// Find if impactor node is inside or outside of the face
  // -------------------------------------------------------

  if(gap < -tolerance)
    is_inside = true;
  else
    is_inside = false;

  // ----------------------------------------------------
  /// Find if impactor node is in projection area of face
  // ----------------------------------------------------

  UInt node_1 = surface_element * nb_nodes_element;
  Real * position_node_1 = &(current_position[connectivity[node_1 + 0] * spatial_dimension]);
  Real * position_node_2 = &(current_position[connectivity[node_1 + 1] * spatial_dimension]);

  Real tmp_vector_1_imp[2];
  Real tmp_vector_1_2[2];

  // find vectors from master node 1 to impactor and master node 2
  Math::vector_2d(position_node_1, position_node_2, tmp_vector_1_2);
  Math::vector_2d(position_node_1, projected_position, tmp_vector_1_imp);

  Real length_difference = Math::norm2(tmp_vector_1_imp) - Math::norm2(tmp_vector_1_2);

  // the projection is outside if the test area is larger than the area of the triangle
  if(length_difference > tolerance)
    is_in_projection_area = false;
  else
    is_in_projection_area = true;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch3dExplicit::checkPenetrationSituationTriangle3(const UInt impactor_node,
								 const UInt surface_element,
								 bool & is_inside,
								 bool & is_in_projection_area) {

  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(spatial_dimension == 3, "wrong spatial dimension (=" << spatial_dimension << ") for checkPenetrationSituationTriangle3");
  const UInt dim = spatial_dimension;
  const ElementType type = _triangle_3;
  const UInt nb_nodes_element = Mesh::getNbNodesPerElement(type);
  const Real tolerance = std::numeric_limits<Real>::epsilon();

  Real * current_position = contact.getModel().getCurrentPosition().values;
  UInt * connectivity = mesh.getConnectivity(type).values;

  Real gap;
  Real normal[3];
  Real projected_position[3];
  computeComponentsOfProjectionTriangle3(impactor_node, surface_element, normal, gap, projected_position);

  // -------------------------------------------------------
  /// Find if impactor node is inside or outside of the face
  // -------------------------------------------------------

  if(gap < -tolerance)
    is_inside = true;
  else
    is_inside = false;

  // ----------------------------------------------------
  /// Find if impactor node is in projection area of face
  // ----------------------------------------------------

  UInt node_1 = surface_element * nb_nodes_element;
  Real * position_node_1 = &(current_position[connectivity[node_1 + 0] * spatial_dimension]);
  Real * position_node_2 = &(current_position[connectivity[node_1 + 1] * spatial_dimension]);
  Real * position_node_3 = &(current_position[connectivity[node_1 + 2] * spatial_dimension]);

  Real triangle_area;
  Real test_area = 0.0;
  Real tmp_vector_1[3];
  Real tmp_vector_2[3];
  Real tmp_vector_3[3];

  // find area of triangle
  Math::vector_3d(position_node_1, position_node_2, tmp_vector_1);
  Math::vector_3d(position_node_1, position_node_3, tmp_vector_2);
  Math::vectorProduct3(tmp_vector_1, tmp_vector_2, tmp_vector_3);
  triangle_area = 0.5 * Math::norm3(tmp_vector_3);

  // compute areas of projected position and two nodes of master triangle
  UInt nb_sub_areas = nb_nodes_element;
  UInt node_order[4];
  node_order[0] = 0; node_order[1] = 1; node_order[2] = 2; node_order[3] = 0;
  for(UInt i=0; i < nb_sub_areas; ++i) {
    position_node_1 = &(current_position[connectivity[node_1 + node_order[i  ]] * spatial_dimension]);
    position_node_2 = &(current_position[connectivity[node_1 + node_order[i+1]] * spatial_dimension]);
    Math::vector_3d(projected_position, position_node_1, tmp_vector_1);
    Math::vector_3d(projected_position, position_node_2, tmp_vector_2);
    Math::vectorProduct3(tmp_vector_1, tmp_vector_2, tmp_vector_3);
    test_area += 0.5 * Math::norm3(tmp_vector_3);
  }

  // the projection is outside if the test area is larger than the area of the triangle
  if((test_area - triangle_area) > tolerance)
    is_in_projection_area = false;
  else
    is_in_projection_area = true;

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
