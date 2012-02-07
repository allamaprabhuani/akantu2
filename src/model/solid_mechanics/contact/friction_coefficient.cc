/**
 * @file   friction_coefficient.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Feb 22 11:08:40 2011
 *
 * @brief  
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
#include "friction_coefficient.hh"
#include "contact_rigid.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
FrictionCoefficient::FrictionCoefficient(ContactRigid & contact,
					 const Surface & master_surface) :
  contact(contact), master_surface(master_surface) {
  AKANTU_DEBUG_IN();

  this->contact.setFrictionCoefficient(this);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
FrictionCoefficient::~FrictionCoefficient() {
  AKANTU_DEBUG_IN();

  this->contact.removeFrictionCoefficient(this->master_surface);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FrictionCoefficient::computeFrictionCoefficient(Vector<Real> & fric_coef) {
  AKANTU_DEBUG_IN();

  UInt nb_fric_coef_elements = fric_coef.getSize();
  Real * fric_coef_val = fric_coef.values;
  
  initializeComputeFricCoef();

  for (UInt i=0; i < nb_fric_coef_elements; ++i)
    fric_coef_val[i] = computeFricCoef(i);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void FrictionCoefficient::setParam(const std::string & key,
				   __attribute__((unused)) const std::string & value) {
  AKANTU_DEBUG_IN();
  
  //if(key == "name") name = std::string(value);
  //else AKANTU_DEBUG_ERROR("Unknown friction coefficient property : " << key);
  AKANTU_DEBUG_ERROR("Unknown friction coefficient property : " << key);

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
