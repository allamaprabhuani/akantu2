/**
 * @file   unique_constant_fric_coef.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Tue Mar 22 13:18:24 2011
 *
 * @brief  implementation of unique constant friction coefficient
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
#include "unique_constant_fric_coef.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
UniqueConstantFricCoef::UniqueConstantFricCoef(ContactRigid & contact, 
					       const Surface & master_surface) :
  FrictionCoefficient(contact, master_surface), friction_coefficient(0.) {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UniqueConstantFricCoef::UniqueConstantFricCoef(ContactRigid & contact, 
					       const Surface & master_surface,
					       const Real friction_coefficient) :
  FrictionCoefficient(contact, master_surface), friction_coefficient(friction_coefficient) {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UniqueConstantFricCoef::~UniqueConstantFricCoef() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void UniqueConstantFricCoef::setParam(const std::string & key, const std::string & value) {
  AKANTU_DEBUG_IN();
  
  std::stringstream sstr(value);
  if(key == "mu") { sstr >> this->friction_coefficient; }
  //else if(key == "E") { sstr >> E; }
  else { FrictionCoefficient::setParam(key, value); }
  
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

