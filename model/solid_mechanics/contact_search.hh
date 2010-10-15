/**
 * @file   contact_search.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 10:43:54 2010
 *
 * @brief  Interface of the search class for contact
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_SEARCH_HH__
#define __AKANTU_CONTACT_SEARCH_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {
  class Contact;
  class ContactNeighborStructure;
}

__BEGIN_AKANTU__

class PenetrationList {
public:
  /// number of penetrating nodes
  UInt nb_nodes;

  /// nodes who have penetrated the master surface
  Vector<UInt> penetrating_nodes;

  /// list of penetrated facets
  Vector<UInt> penetrated_facets_offset;
  Vector<UInt> penetrated_facet;

  /// normal of the penetrated facets
  Vector<Real> facets_normals;
  /// gap between the penetrated facets and the node
  Vector<Real> gap;
  /// position of the node projected on the penetrated facets
  Vector<Real> projected_position;
};

/* -------------------------------------------------------------------------- */


class ContactSearch {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactSearch(const Contact & contact,
		const ContactNeighborStructureType & neighbors_structure_type,
		const ContactSearchID & id = "search_contact");

  virtual ~ContactSearch();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the needed structures
  virtual void initSearch();

  /// build the penetration list
  virtual PenetrationList * findPenetration(const Surface & master_surface) = 0;

  /// update the neighbor structure
  virtual void updateStructure(const Surface & master_surface);

  /// check if the neighbor structure need an update
  virtual bool checkIfUpdateStructureNeeded(const Surface & master_surface);

  /// add a new master surface
  void addMasterSurface(const Surface & master_surface);

  /// remove a master surface
  void removeMasterSurface(const Surface & master_surface);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Contact, contact, const Contact &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// id of the object
  ContactSearchID id;

  /// associated contact class
  const Contact & contact;

  /// type of the neighbors structure to create
  const ContactNeighborStructureType & neighbors_structure_type;

  /// structure used to handle neighbors lists
  std::map<Surface, ContactNeighborStructure *> neighbors_structure;
};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_SEARCH_HH__ */
