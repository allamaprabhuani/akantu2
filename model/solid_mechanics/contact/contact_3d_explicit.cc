/**
 * @file   contact_3d_explicit.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 17:15:08 2010
 *
 * @brief Specialization  of the  contact structure for  3d contact  in explicit
 * time scheme
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
