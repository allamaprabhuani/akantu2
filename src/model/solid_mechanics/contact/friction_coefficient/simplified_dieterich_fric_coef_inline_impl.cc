/**
 * @file   simplified_dieterich_fric_coef_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Tue Mar 22 13:18:24 2011
 *
 * @brief  implementation of simplified dieterich fric coef computation
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
inline Real SimplifiedDieterichFricCoef::computeFricCoef(UInt impactor_node_index) {
  AKANTU_DEBUG_IN();
  
  Real friction_coefficient = 0.;

  const Real tolerance = std::numeric_limits<Real>::epsilon();
  AKANTU_DEBUG_ASSERT(this->v_normalizer > 0., "v_normalizer is negativ ! Cannot compute the velocity term of the friction coefficient");
  AKANTU_DEBUG_ASSERT(this->v_normalizer > tolerance, "v_normalizer is zero ! Cannot compute the velocity term of the friction coefficient");
  AKANTU_DEBUG_ASSERT(this->theta_normalizer > 0., "theta_normalizer is negativ ! Cannot compute the theta term of the friction coefficient");
  AKANTU_DEBUG_ASSERT(this->theta_normalizer > tolerance, "theta_normalizer is zero ! Cannot compute the theta term of the friction coefficient");

  Real * relative_sliding_velocities_val = this->relative_sliding_velocities->values;
  Real velocity_term = this->a_factor * log(relative_sliding_velocities_val[impactor_node_index] / this->v_normalizer + 1);

  Real * theta_state_variables_val = this->theta_state_variables->values;
  Real theta_term = this->b_factor * log(theta_state_variables_val[impactor_node_index] / this->theta_normalizer + 1);

  friction_coefficient = this->mu_zero + velocity_term + theta_term;

  AKANTU_DEBUG_OUT();
  return friction_coefficient;
}
