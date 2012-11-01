/**
 * @file   velocity_weakening_exponential_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Tue Jun 21 15:14:39 2011
 *
 * @brief  implementation of exponential velocity weakening friction coefficient
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
inline Real VelocityWeakeningExponential::computeFricCoef(UInt impactor_node_index) {
  AKANTU_DEBUG_IN();
  
  Real friction_coefficient = 0.;

  Real * relative_sliding_velocities_val;
  if (instant_velocity)
    relative_sliding_velocities_val = this->relative_sliding_velocities->values;
  else
    relative_sliding_velocities_val = this->generalized_sliding_velocities->values;

  if (!instant_velocity && (*node_stick_status)(impactor_node_index))
    friction_coefficient = this->static_friction_coefficient;
  else {
    Real velocity_term = 1 - exp((-1.) * this->alpha * relative_sliding_velocities_val[impactor_node_index]);
  
    friction_coefficient = this->static_friction_coefficient;
    friction_coefficient += (this->dynamic_friction_coefficient - this->static_friction_coefficient) * velocity_term;
  }

  AKANTU_DEBUG_OUT();
  return friction_coefficient;
}

/* -------------------------------------------------------------------------- */
inline void VelocityWeakeningExponential::computeAlpha() {
  AKANTU_DEBUG_IN();
  
  const Real tolerance = std::numeric_limits<Real>::epsilon();
  AKANTU_DEBUG_ASSERT(std::abs(power) > tolerance, "Power is equal zero!! Cannot be the case: division by zero");

  this->alpha = std::sqrt((this->static_friction_coefficient - this->dynamic_friction_coefficient) / this->power);

  AKANTU_DEBUG_OUT();
}
