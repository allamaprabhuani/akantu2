/**
 * @file   contact_neighbor_structure.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 12:36:15 2010
 *
 * @brief  Interface of the structure handling the neighbor lists
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__
#define __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "contact_search.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class NeighborList {
public:
  /// number of impactor nodes
  UInt nb_nodes;

  /// @todo check if data must be stored by element type (different surface elements)

  /// list of nodes of slave surfaces near the master one
  Vector<UInt> impactor_nodes;

  /// neighbor facets (sparse storage)
  Vector<UInt> facets_offset;
  Vector<UInt> facets;
};

/* -------------------------------------------------------------------------- */


class ContactNeighborStructure {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactNeighborStructure(const ContactSearch & contact_search,
			   const Surface & master_surface,
			   const ContactNeighborStructureID & id = "contact_neighbor_structure_id");

  virtual ~ContactNeighborStructure() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the structure
  virtual void init() = 0;
  
  /// update the structure
  virtual void update() = 0;
  
  /// check if an update is needed
  virtual bool check();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the neighbor list
  virtual NeighborList * getNeighborList() = 0;

  AKANTU_GET_MACRO(ContactSearch, contact_search, const ContactSearch &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the object
  ContactNeighborStructureID id;

  /// associated contact search class
  const ContactSearch & contact_search;
  
  /// associated master surface
  const Surface & master_surface;
};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__ */
