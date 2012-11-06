/**
 * @file   ruina_slowness_fric_coef_inline_impl.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Tue Mar 22 13:18:24 2011
 *
 * @brief  implementation of theta formulas
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
template<bool compute_analytic_solution>
inline Real RuinaSlownessFricCoef<compute_analytic_solution>::computeAnalyticTheta(Real previous_theta, Real sliding_speed, Real delta_t) {
  AKANTU_DEBUG_IN();

  Real theta;
  const Real tolerance = 10 * std::numeric_limits<Real>::epsilon();

  Real minus_lambda = sliding_speed / this->d_zero;
  
  if (minus_lambda > tolerance) {
    Real e_term = exp((-minus_lambda) * delta_t);
    theta = (e_term - 1) / (-minus_lambda) + previous_theta * e_term;
  }
  else {
    theta = delta_t + previous_theta;
  }

  AKANTU_DEBUG_OUT();
  return theta;
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
inline Real RuinaSlownessFricCoef<compute_analytic_solution>::computeImplicitTheta(Real previous_theta, Real sliding_speed, Real delta_t) {
  AKANTU_DEBUG_IN();

  Real theta_numerator = delta_t + previous_theta;
  Real theta_denominator = sliding_speed / this->d_zero * delta_t + 1;
  
  Real theta = theta_numerator / theta_denominator;
  
  AKANTU_DEBUG_OUT();
  return theta;
}
