/**
 * @file   contact.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 08 15:20:20 2010
 *
 * @brief  Common part of contact classes
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
#include "contact.hh"
#include "contact_search.hh"
#include "aka_common.hh"
#include "contact_3d_explicit.hh"
#include "contact_search_explicit.hh"
#include "contact_2d_explicit.hh"
#include "contact_search_2d_explicit.hh"
#include "contact_rigid.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Contact::Contact(const SolidMechanicsModel & model,
		 const ContactType & type,
		 const ContactID & id,
		 const MemoryID & memory_id) :
  Memory(memory_id), id(id), model(model),
  type(type),
  node_to_elements_offset("node_to_elements_offset", id, memory_id),
  node_to_elements("node_to_elements", id, memory_id) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Contact::~Contact() {
  AKANTU_DEBUG_IN();

  delete contact_search;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::initContact(bool add_surfaces_flag) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model.getFEEngine().getMesh();

  /// Build surfaces if not done yet
  if(mesh.getNbSurfaces() == 0) { /* initialise nb_surfaces to zero in mesh_io */
    MeshUtils::buildFacets(mesh);
    MeshUtils::buildSurfaceID(mesh);
  }

  /// Detect facet types
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt  nb_facet_types = 0;
  ElementType facet_type[_max_element_type];
  UInt spatial_dimension = mesh.getSpatialDimension();

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(mesh.getSpatialDimension(type) != spatial_dimension) continue;
    facet_type[nb_facet_types] = mesh.getFacetType(type);
    nb_facet_types++;
  }

  CSR<UInt> suface_nodes;
  MeshUtils::buildNodesPerSurface(mesh, suface_nodes);
  suface_nodes.copy(surface_to_nodes_offset, surface_to_nodes);


  for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {
    ElementType type = facet_type[el_type];

    this->node_to_elements_offset.alloc(0, 1, type, _not_ghost);
    this->node_to_elements.alloc(0, 1, type, _not_ghost);

    CSR<UInt> node_to_elem;
    MeshUtils::buildNode2ElementsElementTypeMap(mesh, node_to_elem, type);

    node_to_elem.copy(node_to_elements_offset(type, _not_ghost), node_to_elements(type, _not_ghost));
  }

  if(add_surfaces_flag) {
    for (UInt i = 0; i < mesh.getNbSurfaces(); ++i) {
      addMasterSurface(i);
      std::cout << "Master surface " << i << " has been added automatically" << std::endl;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::initSearch() {
  AKANTU_DEBUG_IN();

  contact_search->initSearch();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::initNeighborStructure() {
  AKANTU_DEBUG_IN();

  contact_search->initNeighborStructure();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::initNeighborStructure(__attribute__ ((unused)) const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  contact_search->initNeighborStructure();

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Contact::checkAndUpdate() {
  AKANTU_DEBUG_IN();

  std::vector<Surface>::iterator it;
  for (it = master_surfaces.begin(); it != master_surfaces.end(); ++it) {
    if(contact_search->checkIfUpdateStructureNeeded(*it)) {
      contact_search->updateStructure(*it);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::updateContact() {
  AKANTU_DEBUG_IN();

  std::vector<Surface>::iterator it;
  for (it = master_surfaces.begin(); it != master_surfaces.end(); ++it) {
    contact_search->updateStructure(*it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::addMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(std::find(master_surfaces.begin(),
				master_surfaces.end(),
				master_surface) == master_surfaces.end(),
		      "Master surface already registered in the master surface list");

  master_surfaces.push_back(master_surface);
  contact_search->addMasterSurface(master_surface);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::removeMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(std::find(master_surfaces.begin(),
				master_surfaces.end(),
				master_surface) != master_surfaces.end(),
		      "Master surface not registered in the master surface list");

  std::remove(master_surfaces.begin(),
	      master_surfaces.end(),
	      master_surface);
  contact_search->removeMasterSurface(master_surface);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Contact * Contact::newContact(const SolidMechanicsModel & model,
			      const ContactType & contact_type,
			      const ContactSearchType & contact_search_type,
			      const ContactNeighborStructureType & contact_neighbor_structure_type,
			      const ContactID & id,
			      const MemoryID & memory_id) {
  AKANTU_DEBUG_IN();

  Contact * tmp_contact = NULL;
  ContactSearch * tmp_search __attribute__ ((unused)) =  NULL;

  switch(contact_type) {
  case _ct_2d_expli: {
    tmp_contact = new Contact2dExplicit(model, contact_type, id, memory_id);
    break;
  }
  case _ct_3d_expli: {
    tmp_contact = new Contact3dExplicit(model, contact_type, id, memory_id);
    break;
  }
  case _ct_rigid: {
    tmp_contact = new ContactRigid(model, contact_type, id, memory_id);
    break;
  }
  case _ct_not_defined: {
    AKANTU_DEBUG_ERROR("Not a valid contact type : " << contact_type);
    break;
  }
  }

  std::stringstream sstr;
  sstr << id << ":contact_search";

  switch(contact_search_type) {
  case _cst_2d_expli: {
    tmp_search = new ContactSearch2dExplicit(*tmp_contact, contact_neighbor_structure_type, contact_search_type, sstr.str());
    break;
  }
  case _cst_expli: {
    tmp_search = new ContactSearchExplicit(*tmp_contact, contact_neighbor_structure_type, contact_search_type, sstr.str());
    break;
  }
  case _cst_not_defined: {
    AKANTU_DEBUG_ERROR("Not a valid contact search type : " << contact_search_type);
    break;
  }
  }

  AKANTU_DEBUG_OUT();
  return tmp_contact;
}

__END_AKANTU__
