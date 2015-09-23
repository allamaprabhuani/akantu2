/**
 * @file   remove_damaged_with damage_rate_weight_function_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief Implementation of inline function of remove damaged with
 * damage rate weight function
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
template<UInt spatial_dimension>
inline void RemoveDamagedWithDamageRateWeightFunction<spatial_dimension>::selectType(__attribute__((unused)) ElementType type1,
										     __attribute__((unused)) GhostType ghost_type1,
										     ElementType type2,
										     GhostType ghost_type2) {
  /// select the damage arrays for given types: For optimization
  selected_damage_with_damage_rate = &(this->material.template getArray<Real>("damage",type2, ghost_type2));
  selected_damage_rate_with_damage_rate = &(this->material.template getArray<Real>("damage-rate",type2, ghost_type2));
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real RemoveDamagedWithDamageRateWeightFunction<spatial_dimension>::operator()(Real r,
										     const __attribute__((unused)) QuadraturePoint & q1,
										     const QuadraturePoint & q2) {
  /// compute the weight
  UInt quad = q2.global_num;

  if(q1.global_num == quad) return 1.;

  Real D = (*selected_damage_with_damage_rate)(quad);
  Real w = 0.;
  Real alphaexp = 1.;
  Real betaexp = 2.;
  if(D < damage_limit_with_damage_rate) {
    Real alpha = std::max(0., 1. - pow((r*r / this->R2),alphaexp));
    w = pow(alpha, betaexp);
  }

  return w;
}

