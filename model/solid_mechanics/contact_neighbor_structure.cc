/**
 * @file   contact_neighbor_structure.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 12:42:26 2010
 *
 * @brief  
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
#include "contact_neighbor_structure.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NeighborList::NeighborList() : impactor_nodes(Vector<UInt>(0, 1, "impactors")) {
  AKANTU_DEBUG_IN();

  for (UInt i = 0; i < _max_element_type; ++i) {
    facets_offset[i] = NULL;
    facets       [i] = NULL;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NeighborList::~NeighborList() {
  AKANTU_DEBUG_IN();

  for (UInt i = 0; i < _max_element_type; ++i) {
    if(facets_offset[i]) delete facets_offset[i];
    if(facets       [i]) delete facets       [i];
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactNeighborStructure::ContactNeighborStructure(const ContactSearch & contact_search,
						   const Surface & master_surface,
						   const ContactNeighborStructureType & type,
						   const ContactNeighborStructureID & id) :
  id(id), contact_search(contact_search), master_surface(master_surface),
  type(type) {
  AKANTU_DEBUG_IN();

  neighbor_list = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactNeighborStructure::~ContactNeighborStructure() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool ContactNeighborStructure::check() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ERROR("Check not implemented for this neighbors structure : " << id);
  AKANTU_DEBUG_OUT();

  return false;
}

__END_AKANTU__
