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
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact.hh"
#include "contact_search.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Contact::Contact(const SolidMechanicsModel & model,
		 ContactSearch & contact_search,
		 const ContactID & id,
		 const MemoryID & memory_id) :
  Memory(memory_id), id(id), model(model), contact_search(&contact_search) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

Contact::~Contact() {
  AKANTU_DEBUG_IN();

  delete contact_search;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Contact::initContact() {
  AKANTU_DEBUG_IN();
  
  contact_search->initSearch();

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
			      const ContactID & id) {
  AKANTU_DEBUG_IN();
  
  std::stringstream sstr;
  sstr << id << ":contact_search";

  ContactSearch * tmp_search = NULL;
  switch(contact_search_type) {
  case _cst_not_defined:
    //    tmp_search = new ContactSearch(this, contact_neighbor_structure_type, sstr.str());
    AKANTU_DEBUG_ERROR("Not a valid contact search type : " << contact_search_type);
    break;
  }

  Contact * tmp_contact = NULL;
  switch(contact_search_type) {
  case _ct_not_defined:
    // tmp_contact = new Contact(model, tmp_search, id);
    AKANTU_DEBUG_ERROR("Not a valid contact type : " << contact_type);
    break;
  }

  AKANTU_DEBUG_OUT();
  return tmp_contact;
}

__END_AKANTU__
