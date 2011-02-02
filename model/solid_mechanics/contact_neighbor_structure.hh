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
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

#ifndef __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__
#define __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "mesh.hh"
#include "contact_search.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class NeighborList {
public:
  NeighborList();
  virtual ~NeighborList();
public:
  /// number of impactor nodes
  //UInt nb_nodes;

  /// list of nodes of slave surfaces near the master one
  Vector<UInt> impactor_nodes;

  /// neighbor facets (sparse storage)
  ByElementTypeUInt facets_offset;
  ByElementTypeUInt facets;
};

/* -------------------------------------------------------------------------- */


class ContactNeighborStructure {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ContactNeighborStructure(const ContactSearch & contact_search,
			   const Surface & master_surface,
			   const ContactNeighborStructureType & type,
			   const ContactNeighborStructureID & id = "contact_neighbor_structure_id");

  virtual ~ContactNeighborStructure();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the structure
  virtual void initNeighborStructure() = 0;
  
  /// update the structure
  virtual void update() = 0;
  
  /// check if an update is needed
  virtual bool check();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the neighbor list
  AKANTU_GET_MACRO(NeighborList, *neighbor_list, const NeighborList &);

  AKANTU_GET_MACRO(ID, id, const ContactNeighborStructureID)

  AKANTU_GET_MACRO(ContactSearch, contact_search, const ContactSearch &);

  AKANTU_GET_MACRO(MasterSurface, master_surface, const Surface &);

  AKANTU_GET_MACRO(Type, type, const ContactNeighborStructureType &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the object
  ContactNeighborStructureID id;

  /// associated contact search class
  const ContactSearch & contact_search;
  
  /// associated master surface
  Surface master_surface;

  /// neighbor list
  NeighborList * neighbor_list;

  /// type of contact neighbor structure object
  ContactNeighborStructureType type;

};

__END_AKANTU__

#endif /* __AKANTU_CONTACT_NEIGHBOR_STRUCTURE_HH__ */
