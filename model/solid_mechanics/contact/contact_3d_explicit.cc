/**
 * @file   contact_3d_explicit.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 17:15:08 2010
 *
 * @brief  Specialization of the contact structure for 3d contact in explicit time scheme
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_3d_explicit.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
Contact3dExplicit::Contact3dExplicit(const SolidMechanicsModel & model,
				     const ContactType & type,
				     const ContactID & id,
				     const MemoryID & memory_id) :
  Contact(model, type, id, memory_id) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

Contact3dExplicit::~Contact3dExplicit() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void Contact3dExplicit::solveContact() {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */





__END_AKANTU__
