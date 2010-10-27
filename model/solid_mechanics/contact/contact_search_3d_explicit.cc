/**
 * @file   contact_search_3d_explicit.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 18:49:04 2010
 *
 * @brief  Specialization of the contact search structure for 3D within an explicit time scheme 
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_search_3d_explicit.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactSearch3dExplicit::ContactSearch3dExplicit(Contact & contact,
						 const ContactNeighborStructureType & neighbors_structure_type,
						 const ContactSearchID & id) :
  ContactSearch(contact, neighbors_structure_type, id) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
PenetrationList * ContactSearch3dExplicit::findPenetration(const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
  return ;
}

/* -------------------------------------------------------------------------- */




__END_AKANTU__
