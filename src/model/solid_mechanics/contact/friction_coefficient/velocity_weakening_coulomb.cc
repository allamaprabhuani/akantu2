/**
 * @file   velocity_weakening_coulomb.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Wed May 11 15:14:03 2011
 *
 * @brief  implementation of velocity weakening constant friction coefficient
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
#include "velocity_weakening_coulomb.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
VelocityWeakeningCoulomb::VelocityWeakeningCoulomb(ContactRigid & contact, 
						   const Surface & master_surface) :
  FrictionCoefficient(contact, master_surface), static_friction_coefficient(0.), 
  dynamic_friction_coefficient(0.) {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
VelocityWeakeningCoulomb::VelocityWeakeningCoulomb(ContactRigid & contact, 
						   const Surface & master_surface,
						   const Real static_friction_coefficient,
						   const Real dynamic_friction_coefficient) :
  FrictionCoefficient(contact, master_surface), 
  static_friction_coefficient(static_friction_coefficient),
  dynamic_friction_coefficient(dynamic_friction_coefficient){
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
VelocityWeakeningCoulomb::~VelocityWeakeningCoulomb() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void VelocityWeakeningCoulomb::setParam(const std::string & key, const std::string & value) {
  AKANTU_DEBUG_IN();
  
  std::stringstream sstr(value);
  if(key == "mu_s") { sstr >> this->static_friction_coefficient; }
  else  if(key == "mu_d") { sstr >> this->dynamic_friction_coefficient; }
  else { FrictionCoefficient::setParam(key, value); }
  
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

