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
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_neighbor_structure.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactNeighborStructure::ContactNeighborStructure(const ContactSearch & contact_search,
						   const Surface & master_surface,
						   const ContactNeighborStructureID & id) :
  id(id), contact_search(contact_search), master_surface(master_surface) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

bool ContactNeighborStructure::check() {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ERROR("Check not implemented for this neighbors structure : " << id);
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
