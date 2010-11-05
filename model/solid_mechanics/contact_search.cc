/**
 * @file   contact_search.cc
 * @author David Kammer <kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 11:46:34 2010
 *
 * @brief
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_search.hh"
#include "contact.hh"
#include "contact_neighbor_structure.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactSearch::ContactSearch(Contact & contact,
			     const ContactNeighborStructureType & neighbors_structure_type,
			     const ContactSearchID & id) :
  id(id), contact(contact), neighbors_structure_type(neighbors_structure_type) {
  AKANTU_DEBUG_IN();

  contact.setContactSearch(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactSearch::~ContactSearch() {
  AKANTU_DEBUG_IN();

  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    delete it->second;
  }

  neighbors_structure.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::initSearch() {
  AKANTU_DEBUG_IN();

  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    it->second->init();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::addMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) == neighbors_structure.end(),
		      "Master surface already registered in the search object " << id);

  ContactNeighborStructure * tmp_neighbors_structure = NULL;

  std::stringstream sstr;
  sstr << id << ":contact_neighbor_structure:" << neighbors_structure_type << ":" << master_surface;

  switch(neighbors_structure_type) {
  case _cnst_regular_grid :
    // if mesh.getSpatialDimension() == 2 then RegularGridNeighborStructure<2>(...);
    // else if mesh.getSpatialDimension() == 3 then RegularGridNeighborStructure<3>(...);
    // else error
    break;
  case _cnst_not_defined :
    //    tmp_neighbors_structure = new ContactNeighborStructureGrid2d(this, master_surface, sstr.str());
    AKANTU_DEBUG_ERROR("Not a valid neighbors structure type : " << neighbors_structure_type);
    break;
  }

  neighbors_structure[master_surface] = tmp_neighbors_structure;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::removeMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  delete neighbors_structure[master_surface];
  neighbors_structure.erase(neighbors_structure.find(master_surface));

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch::updateStructure(const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  neighbors_structure[master_surface]->update();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool ContactSearch::checkIfUpdateStructureNeeded(const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  bool check = neighbors_structure[master_surface]->check();

  AKANTU_DEBUG_OUT();
  return check;
}



__END_AKANTU__
