/**
 * @file   unique_constant_fric_coef_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Tue Mar 22 13:18:24 2011
 *
 * @brief  implementation of a unique constant friction coefficient
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
inline Real UniqueConstantFricCoef::computeFricCoef(__attribute__((unused)) UInt impactor_node_index) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
  return this->friction_coefficient;
}
