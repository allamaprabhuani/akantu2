/**
 * @file   velocity_weakening_exponential.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Jun 20 15:02:47 2011
 *
 * @brief  implementation of an exponential velocity weakening friction coef
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
#include "velocity_weakening_exponential.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
VelocityWeakeningExponential::VelocityWeakeningExponential(ContactRigid & contact, 
							   const Surface & master_surface,
							   const Real static_friction_coefficient,
							   const Real dynamic_friction_coefficient,
							   const Real power) :
  FrictionCoefficient(contact, master_surface),
  VelocityDependentFricCoef(contact, master_surface), 
  HistoricVelocityFricCoef(contact, master_surface),
  static_friction_coefficient(static_friction_coefficient),
  dynamic_friction_coefficient(dynamic_friction_coefficient),
  power(power),
  instant_velocity(true) {
  AKANTU_DEBUG_IN();

  computeAlpha();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
VelocityWeakeningExponential::VelocityWeakeningExponential(ContactRigid & contact, 
							   const Surface & master_surface,
							   const Real static_friction_coefficient,
							   const Real dynamic_friction_coefficient,
							   const Real power,
							   const Real beta) :
  FrictionCoefficient(contact, master_surface),
  VelocityDependentFricCoef(contact, master_surface), 
  HistoricVelocityFricCoef(contact, master_surface, beta), 
  static_friction_coefficient(static_friction_coefficient),
  dynamic_friction_coefficient(dynamic_friction_coefficient),
  power(power),
  instant_velocity(false) {
  AKANTU_DEBUG_IN();

  computeAlpha();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
VelocityWeakeningExponential::~VelocityWeakeningExponential() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void VelocityWeakeningExponential::initializeComputeFricCoef() {
  AKANTU_DEBUG_IN();
  
  // compute relative sliding velocities
  if (instant_velocity)
    VelocityDependentFricCoef::initializeComputeFricCoef();
  else
    HistoricVelocityFricCoef::initializeComputeFricCoef();

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__

