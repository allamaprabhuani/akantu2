/**
 * @file   newmark-beta_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep 30 11:45:20 2010
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
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
/*
 * @f$ \tilde{u_{n+1}} = u_{n} +  \Delta t \dot{u}_n + \frac{\Delta t^2}{2} \ddot{u}_n @f$
 * @f$ \tilde{\dot{u}_{n+1}} = \dot{u}_{n} +  \Delta t \ddot{u}_{n} @f$
 * @f$ \tilde{\ddot{u}_{n}} = \ddot{u}_{n} @f$
 */
void NewmarkBeta::integrationSchemePred(Real delta_t,
					Vector<Real> & u,
					Vector<Real> & u_dot,
					Vector<Real> & u_dot_dot,
					Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      Real dt_a_n = delta_t * *u_dot_dot_val;

      *u_val += delta_t * *u_dot_val + delta_t / 2. * dt_a_n;
      *u_dot_val += dt_a_n;
      // *u_dot_dot_val = *u_dot_dot_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
 * @f$ u_{n+1} = \tilde{u_{n+1}} + \alpha \beta \Delta t^2 \delta \ddot{u}_n @f$
 * @f$ \dot{u}_{n+1} = \tilde{\dot{u}_{n+1}} + \beta \Delta t * \delta \ddot{u}_{n+1} @f$
 * @f$ \ddot{u}_{n+1} = \tilde{\ddot{u}_{n+1}} + \delta \ddot{u}_{n+1} @f$
 */
void NewmarkBeta::integrationSchemeCorrAccel(Real delta_t,
					    Vector<Real> & u,
					    Vector<Real> & u_dot,
					    Vector<Real> & u_dot_dot,
					    Vector<bool> & boundary,
					    Vector<Real> & delta
					    ) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<_acceleration_corrector>(delta_t,
						 u,
						 u_dot,
						 u_dot_dot,
						 boundary,
						 delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
 * @f$ u_{n+1} = \tilde{u_{n+1}} +  \alpha \Delta t \delta \dot{u}_n @f$
 * @f$ \dot{u}_{n+1} = \tilde{\dot{u}_{n+1}} + \delta \dot{u}_{n+1} @f$
 * @f$ \ddot{u}_{n+1} = \tilde{\ddot{u}_{n+1}} + \frac{1}{\beta \Delta t} \delta \dot{u}_{n+1} @f$
 */
void NewmarkBeta::integrationSchemeCorrVeloc(Real delta_t,
					     Vector<Real> & u,
					     Vector<Real> & u_dot,
					     Vector<Real> & u_dot_dot,
					     Vector<bool> & boundary,
					     Vector<Real> & delta
					     ) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<_velocity_corrector>(delta_t,
					     u,
					     u_dot,
					     u_dot_dot,
					     boundary,
					     delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
 * @f$ u_{n+1} = \tilde{u_{n+1}} + \delta u_n @f$
 * @f$ \dot{u}_{n+1} = \tilde{\dot{u}_{n+1}} +  \frac{1}{\alpha \Delta t} \delta u_{n+1} @f$
 * @f$ \ddot{u}_{n+1} = \tilde{\ddot{u}_{n+1}} + \frac{1}{\alpha \beta \Delta t^2} \delta u_{n+1} @f$
 */
void NewmarkBeta::integrationSchemeCorrDispl(Real delta_t,
					     Vector<Real> & u,
					     Vector<Real> & u_dot,
					     Vector<Real> & u_dot_dot,
					     Vector<bool> & boundary,
					     Vector<Real> & delta
					     ) {
  AKANTU_DEBUG_IN();

  integrationSchemeCorr<_displacement_corrector>(delta_t,
						 u,
						 u_dot,
						 u_dot_dot,
						 boundary,
						 delta);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<>
Real NewmarkBeta::getAccelerationCoefficient<NewmarkBeta::_acceleration_corrector>(Real delta_t) {
  return 1.;
}
template<>
Real NewmarkBeta::getAccelerationCoefficient<NewmarkBeta::_velocity_corrector>(Real delta_t) {
  return 1. / (beta * delta_t);
}
template<>
Real NewmarkBeta::getAccelerationCoefficient<NewmarkBeta::_displacement_corrector>(Real delta_t) {
  return 1. / (alpha * beta * delta_t * delta_t);
}

/* -------------------------------------------------------------------------- */
template<>
Real NewmarkBeta::getVelocityCoefficient<NewmarkBeta::_acceleration_corrector>(Real delta_t) {
  return beta * delta_t;
}
template<>
Real NewmarkBeta::getVelocityCoefficient<NewmarkBeta::_velocity_corrector>(Real delta_t) {
  return 1.;
}
template<>
Real NewmarkBeta::getVelocityCoefficient<NewmarkBeta::_displacement_corrector>(Real delta_t) {
  return 1. / (alpha * delta_t);
}

/* -------------------------------------------------------------------------- */
template<>
Real NewmarkBeta::getDisplacementCoefficient<NewmarkBeta::_acceleration_corrector>(Real delta_t) {
  return alpha * beta * delta_t * delta_t;
}
template<>
Real NewmarkBeta::getDisplacementCoefficient<NewmarkBeta::_velocity_corrector>(Real delta_t) {
  return alpha * delta_t;
}
template<>
Real NewmarkBeta::getDisplacementCoefficient<NewmarkBeta::_displacement_corrector>(Real delta_t) {
  return 1.;
}




/* -------------------------------------------------------------------------- */
template<NewmarkBeta::IntegrationSchemeCorrectorType type>
void NewmarkBeta::integrationSchemeCorr(Real delta_t,
					Vector<Real> & u,
					Vector<Real> & u_dot,
					Vector<Real> & u_dot_dot,
					Vector<bool> & boundary,
					Vector<Real> & delta
					) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real c = getAccelerationCoefficient<type>(delta_t);
  Real d = getVelocityCoefficient<type>(delta_t);
  Real e = getDisplacementCoefficient<type>(delta_t);

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  Real * delta_val     = delta.values;
  bool * boundary_val  = boundary.values;

  for (UInt dof = 0; dof < nb_degre_of_freedom; dof++) {
    if(!(*boundary_val)) {
      *u_val         += e * *delta_val;
      *u_dot_val     += d * *delta_val;
      *u_dot_dot_val += c * *delta_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    delta_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}
