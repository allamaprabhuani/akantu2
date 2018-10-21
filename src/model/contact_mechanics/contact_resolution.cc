/**
 * @file  contact_resolution.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Sep 20 2018
 * @date last modification: Thu Sep 20 2018
 *
 * @brief  Implementation of the base class Contact ContactResolution
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_resolution.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactResolution::ContactResolution(const ID & id) {

  
}

/* -------------------------------------------------------------------------- */
ContactResolution::~ContactResolution() = default;

/* -------------------------------------------------------------------------- */
/*void ContactResolution::checkIfTypeIsSupported() {
  if (this->supported_type.find(this->contact_resolution_type) ==
      this->supported_type.end()) {
    AKANTU_EXCEPTION("The contact resolution method"
		     << this->non_linear_solver_type
		     << " is not implemented in the contact resolution "
		     << this->id << "!");
  }
}*/

/* -------------------------------------------------------------------------- */
void ContactResolution::resolutionStep() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();

}
  
}

