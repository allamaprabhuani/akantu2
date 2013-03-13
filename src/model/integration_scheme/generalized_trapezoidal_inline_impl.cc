/**
 * @file   generalized_trapezoidal_inline_impl.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Jul 04 15:07:57 2011
 *
 * @brief  implementation of inline functions
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

void GeneralizedTrapezoidal::integrationSchemePred(Real delta_t,
						   Array<Real> & u,
						   Array<Real> & u_dot,
						   Array<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degree_of_freedom; d++) {
    if(!(*boundary_val)) {
      *u_val += delta_t * *u_dot_val;
    }
    u_val++;
    u_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::integrationSchemeCorrTemp(Real delta_t,
						       Array<Real> & u,
						       Array<Real> & u_dot,
						       Array<bool> & boundary,
						       Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<GeneralizedTrapezoidal::_temperature_corrector>(delta_t,
									u,
									u_dot,
									boundary,
									delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GeneralizedTrapezoidal::integrationSchemeCorrTempRate(Real delta_t,
							   Array<Real> & u,
							   Array<Real> & u_dot,
							   Array<bool> & boundary,
							   Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<GeneralizedTrapezoidal::_temperature_rate_corrector>(delta_t,
									     u,
									     u_dot,
									     boundary,
									     delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
Real GeneralizedTrapezoidal::getTemperatureCoefficient<GeneralizedTrapezoidal::_temperature_corrector>(__attribute__ ((unused)) Real delta_t) {
  return 1.;
}

template<>
Real GeneralizedTrapezoidal::getTemperatureRateCoefficient<GeneralizedTrapezoidal::_temperature_corrector>(Real delta_t) {
  return 1./(alpha * delta_t);
}

/* -------------------------------------------------------------------------- */
template<>
Real GeneralizedTrapezoidal::getTemperatureCoefficient<GeneralizedTrapezoidal::_temperature_rate_corrector>(Real delta_t) {
  return alpha * delta_t;
}

template<>
Real GeneralizedTrapezoidal::getTemperatureRateCoefficient<GeneralizedTrapezoidal::_temperature_rate_corrector>(__attribute__ ((unused)) Real delta_t) {
  return 1.;
}

/* -------------------------------------------------------------------------- */
template<GeneralizedTrapezoidal::IntegrationSchemeCorrectorType type>
void GeneralizedTrapezoidal::integrationSchemeCorr(Real delta_t,
						   Array<Real> & u,
						   Array<Real> & u_dot,
						   Array<bool> & boundary,
						   Array<Real> & delta) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degree_of_freedom = u.getNbComponent() * nb_nodes;

  Real e = getTemperatureCoefficient<type>(delta_t);
  Real d = getTemperatureRateCoefficient<type>(delta_t);

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * delta_val     = delta.values;
  bool * boundary_val  = boundary.values;

  for (UInt dof = 0; dof < nb_degree_of_freedom; dof++) {
    if(!(*boundary_val)) {
      *u_val         += e * *delta_val;
      *u_dot_val     += d * *delta_val;
    }
    u_val++;
    u_dot_val++;
    delta_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}
