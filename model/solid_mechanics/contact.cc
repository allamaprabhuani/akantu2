/**
 * @file   contact.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 14:55:42 2010
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
#include "contact_search_3d_explicit.hh"
#include "contact_2d_explicit.hh"
#include "contact_search_2d_explicit.hh"
#include "contact_rigid_no_friction.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Contact::Contact(const SolidMechanicsModel & model,
		 const ContactType & type,
		 const ContactID & id,
		 const MemoryID & memory_id) :
  Memory(memory_id), id(id), model(model), friction_coefficient(0.),
  master_surfaces(NULL), surface_to_nodes(NULL), surface_to_nodes_offset(NULL),
  type(type) {

  AKANTU_DEBUG_IN();

  for (UInt i = 0; i < _max_element_type; ++i) {
    node_to_elements_offset[i] = NULL;
    node_to_elements[i] = NULL;
  }

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

  Mesh & mesh = model.getFEM().getMesh();
  const UInt nb_nodes = mesh.getNbNodes();

  /// Build surfaces if not done yet
  if(mesh.getNbSurfaces() == 0) { /* initialise nb_surfaces to zero in mesh_io */
    MeshUtils::buildFacets(mesh,1,0);
    MeshUtils::buildSurfaceID(mesh);
  }

  /// Detect facet types
  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  UInt nb_types = type_list.size();
  UInt  nb_facet_types = 0;
  ElementType facet_type[_max_element_type];
  UInt spatial_dimension = mesh.getSpatialDimension();

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(mesh.getSpatialDimension(type) != spatial_dimension) continue;
    facet_type[nb_facet_types] = mesh.getFacetElementType(type);
    nb_facet_types++;
  }

  const UInt nb_surfaces = mesh.getNbSurfaces();

  UInt * surf_nodes_id = new UInt [nb_nodes*nb_surfaces];
  memset(surf_nodes_id, 0, nb_nodes*nb_surfaces*sizeof(Int));

  /// Fill class members node_to_elements_offset and node_to_element
  for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {

    ElementType type = facet_type[el_type];
    const UInt *surf_id_val = mesh.getSurfaceId(type).values;
    std::stringstream sstr_name_1; sstr_name_1 << id << ":node_to_elements_offset:" << type;
    this->node_to_elements_offset[type] = &(alloc<UInt>(sstr_name_1.str(),0,1));
    std::stringstream sstr_name_2; sstr_name_2 << id << ":node_to_elements:" << type;
    this->node_to_elements[type] = &(alloc<UInt>(sstr_name_2.str(),0,1));
    MeshUtils::buildNode2ElementsByElementType(mesh, type, *(node_to_elements_offset[type]), *(node_to_elements[type]));
    UInt * node_off_val = node_to_elements_offset[type]->values;
    UInt * node_elem_val = node_to_elements[type]->values;
    for (UInt n = 0; n < nb_nodes; ++n)
      if(node_off_val[n] != node_off_val[n+1])
	for (UInt c_el = node_off_val[n]; c_el < node_off_val[n+1]; ++c_el) {
	  UInt surf_id = surf_id_val[node_elem_val[c_el]];
	  surf_nodes_id[surf_id*nb_nodes + n] = 1;
	}
  }

  /// Count nodes per surfaces
  this->surface_to_nodes_offset.resize(nb_surfaces+1);
  UInt * surf_nodes_off_val = surface_to_nodes_offset.values;
  memset(surf_nodes_off_val, 0, (nb_surfaces + 1)*sizeof(UInt));

  for (UInt i = 0; i < nb_surfaces; ++i)
    for (UInt n = 0; n < nb_nodes; ++n)
      if(surf_nodes_id[i*nb_nodes+n] == 1)
	surf_nodes_off_val[i]++;

  /// convert the occurrence array in a csr one
  for (UInt i = 1; i < nb_surfaces; ++i) surf_nodes_off_val[i] += surf_nodes_off_val[i-1];
  for (UInt i = nb_surfaces; i > 0; --i) surf_nodes_off_val[i] = surf_nodes_off_val[i-1];
  surf_nodes_off_val[0] = 0;

  // std::stringstream stn_name; stn_name << id << ":surface_to_nodes";
  // this->surface_to_nodes = alloc<UInt>(stn_name.str(),surf_nodes_off_val[nb_surfaces],1);

  /// Fill surface_to_nodes (save node index)
  surface_to_nodes.resize(surf_nodes_off_val[nb_surfaces]);
  UInt * surf_nodes_val = surface_to_nodes.values;

  for (UInt i = 0; i < nb_surfaces; ++i)
    for (UInt n = 0; n < nb_nodes; ++n)
      if(surf_nodes_id[i*nb_nodes+n] == 1) {
	surf_nodes_val[surf_nodes_off_val[i]] = n;
	surf_nodes_off_val[i]++;
      }

  /// rearrange surface_to_nodes_offset to start with 0
  for (UInt i = nb_surfaces; i > 0; --i) surf_nodes_off_val[i] = surf_nodes_off_val[i-1];
  surf_nodes_off_val[0] = 0;

  delete [] surf_nodes_id;

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
void Contact::initNeighborStructure(const Surface & master_surface) {
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
  ContactSearch * tmp_search = NULL;

  switch(contact_type) {
  case _ct_2d_expli: { 
    tmp_contact = new Contact2dExplicit(model, contact_type, id, memory_id); 
    break;
  }
  case _ct_3d_expli: {
    tmp_contact = new Contact3dExplicit(model, contact_type, id, memory_id);
    break;
  }
  case _ct_rigid_no_fric: {
    tmp_contact = new ContactRigidNoFriction(model, contact_type, id, memory_id);
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
  case _cst_3d_expli: {
    tmp_search = new ContactSearch3dExplicit(*tmp_contact, contact_neighbor_structure_type, contact_search_type, sstr.str());
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
