/**
 * @file   contact_search.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Oct 08 15:20:20 2010
 *
 * @brief  Interface of the search class for contact
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

#ifndef __AKANTU_CONTACT_SEARCH_HH__
#define __AKANTU_CONTACT_SEARCH_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "contact.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
  class ContactNeighborStructure;
}

__BEGIN_AKANTU__

class PenetrationList {
public:
  PenetrationList(const ID & id);
  virtual ~PenetrationList();
public:
  ID id;

  /// nodes who have penetrated the master surface
  Array<UInt> penetrating_nodes;

  /// list of penetrated facets
  ByElementTypeUInt penetrated_facets_offset;
  ByElementTypeUInt penetrated_facets;

  /// normal of the penetrated facets
  ByElementTypeReal facets_normals;
  /// gap between the penetrated facets and the node
  ByElementTypeReal gaps;
  /// position of the node projected on the penetrated facets
  ByElementTypeReal projected_positions;    
};

/* -------------------------------------------------------------------------- */


class ContactSearch : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ContactSearch(Contact & contact,
		const ContactNeighborStructureType & neighbors_structure_type,
		const ContactSearchType & type,
		const ContactSearchID & id = "search_contact");

  virtual ~ContactSearch();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the needed structures
  virtual void initSearch();

  /// initialize all neighbor structures
  virtual void initNeighborStructure();

  /// initialize one neighbor structure
  virtual void initNeighborStructure(const Surface & master_surface);

  /// build the penetration list
  virtual void findPenetration(const Surface & master_surface, PenetrationList & penetration_list) = 0;

  /// update the neighbor structure
  virtual void updateStructure(const Surface & master_surface);

  /// check if the neighbor structure need an update
  virtual bool checkIfUpdateStructureNeeded(const Surface & master_surface);

  /// add a new master surface
  void addMasterSurface(const Surface & master_surface);

  /// remove a master surface
  void removeMasterSurface(const Surface & master_surface);

private:
  /// compute the maximal increment for all surface nodes in each direction
  void computeMaxIncrement(Real * max_increment);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(ID, id, const ContactSearchID &);

  AKANTU_GET_MACRO(Contact, contact, const Contact &);

  AKANTU_GET_MACRO(Type, type, const ContactSearchType &);

  const ContactNeighborStructure & getContactNeighborStructure(const Surface & master_surface) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the object
  ContactSearchID id;

  /// associated contact class
  const Contact & contact;

  /// type of the neighbors structure to create
  const ContactNeighborStructureType & neighbors_structure_type;

  /// structure used to handle neighbors lists
  std::map<Surface, ContactNeighborStructure *> neighbors_structure;

  /// type of contact search object
  ContactSearchType type;
};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_SEARCH_HH__ */
