/**
 * @file   damaged_weight_function_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Fri Apr 13 2012
 * @date last modification: Thu Jun 05 2014
 *
 * @brief Implementation of inline function of damaged weight function 
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
inline Real DamagedWeightFunction::operator()(Real r,
					      const __attribute__((unused)) QuadraturePoint & q1,
					      const QuadraturePoint & q2) {
  /// compute the weight
  UInt quad = q2.global_num;
  Array<Real> & dam_array = (*this->damage)(q2.type, q2.ghost_type);
  Real D = dam_array(quad);
  Real Radius_t = 0;
  Real Radius_init = this->R2;

  //    if(D <= 0.5)
  //      {
  //	Radius_t = 2*D*Radius_init;
  //      }
  //    else
  //      {
  //	Radius_t = 2*Radius_init*(1-D);
  //      }
  //

  Radius_t = Radius_init*(1-D);


  Radius_init *= Radius_init;
  Radius_t *= Radius_t;

  if(Radius_t < Math::getTolerance()) {
    Radius_t = 0.001*Radius_init;
  }

  Real expb = (2*std::log(0.51))/(std::log(1.0-0.49*Radius_t/Radius_init));
  Int  expb_floor=std::floor(expb);
  Real b = expb_floor + expb_floor%2;
  Real alpha = std::max(0., 1. - r*r / Radius_init);
  Real w = std::pow(alpha,b);
  return w;
}

/* -------------------------------------------------------------------------- */
inline void DamagedWeightFunction::init() {
  this->damage = &(this->manager.registerWeightFunctionInternal("damage"));
}
