/**
 * @file   rice_kuwano_modified.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Fri Nov 04 14:28:28 2011
 *
 * @brief  implementation of the functions
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
#include "rice_kuwano_modified.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
RiceKuwanoModified::RiceKuwanoModified(ContactRigid & contact, 
				       const Surface & master_surface,
				       const Real static_friction_coefficient,
				       const Real dynamic_friction_coefficient,
				       const Real reference_velocity,
				       const Real alpha) :
  FrictionCoefficient(contact, master_surface),
  VelocityDependentFricCoef(contact, master_surface), 
  static_friction_coefficient(static_friction_coefficient),
  dynamic_friction_coefficient(dynamic_friction_coefficient),
  reference_velocity(reference_velocity),
  alpha(alpha) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
RiceKuwanoModified::~RiceKuwanoModified() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void RiceKuwanoModified::initializeComputeFricCoef() {
  AKANTU_DEBUG_IN();
  
  // compute relative sliding velocities
  VelocityDependentFricCoef::initializeComputeFricCoef();

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__

