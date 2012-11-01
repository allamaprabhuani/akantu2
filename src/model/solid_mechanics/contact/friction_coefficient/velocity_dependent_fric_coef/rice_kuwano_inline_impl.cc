/**
 * @file   rice_kuwano_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Fri Nov 11 16:56:08 2011
 *
 * @brief  implementation of inlined functions
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
inline Real RiceKuwano::computeFricCoef(UInt impactor_node_index) {
  AKANTU_DEBUG_IN();
  
  Real friction_coefficient = 0.;

  Real rel_sliding_vel = (*relative_sliding_velocities)(impactor_node_index);

  if (rel_sliding_vel < this->reference_velocity)
    friction_coefficient = this->static_friction_coefficient;
  else {
    friction_coefficient = (this->static_friction_coefficient - this->dynamic_friction_coefficient) * reference_velocity / rel_sliding_vel;
    friction_coefficient += this->dynamic_friction_coefficient;
    friction_coefficient += this->alpha * rel_sliding_vel;
  }

  AKANTU_DEBUG_OUT();
  return friction_coefficient;
}
