/**
 * @file   newmark-beta_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep 30 11:45:20 2010
 *
 * @brief  implementation of the newmark-@f$\beta@f$ integration scheme
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
void NewmarkBeta::integrationSchemePred(Real delta_t,
					Vector<Real> & u,
					Vector<Real> & u_dot,
					Vector<Real> & u_dot_dot,
					Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real delta_t_2_d2 = delta_t * delta_t / 2.;
  //  Real delta_t_d2   = delta_t / 2.;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      /// @f$ \tilde{u_{n+1}} = u_{n} +  \Delta t \dot{u}_n + (1 - 2 \beta) \frac{\Delta t^2}{2} \ddot{u}_n @f$
      *u_val += delta_t * *u_dot_val + delta_t_2_d2 * (1 - 2 * beta)* *u_dot_dot_val;

      /// @f$ \tilde{\dot{u}_{n+1}} = \dot{u}_{n} +  (1 - \gamma) \Delta t \ddot{u}_{n} @f$
      *u_dot_val += (1 - gamma) * delta_t * *u_dot_dot_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NewmarkBeta::integrationSchemeCorr(Real delta_t,
					Vector<Real> & u,
					Vector<Real> & u_dot,
					Vector<Real> & u_dot_dot,
					Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  //  Real delta_t_d2 =  delta_t / 2.;
  Real delta_t_2_d2 = delta_t * delta_t / 2.;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      /// @f$ u_{n+1} = \tilde{u_{n+1}} + 2 \beta \frac{\Delta t^2}{2} \ddot{u}_n @f$
      *u_val += 2 * beta * delta_t_2_d2 * *u_dot_dot_val;

      /// @f$ \dot{u}_{n+1} = \tilde{\dot{u}_{n+1}} + \gamma \frac{\Delta t}{2} * \ddot{u}_{n+1} @f$
      *u_dot_val += gamma * delta_t * *u_dot_dot_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NewmarkBeta::integrationSchemePredImplicit(Real delta_t,
						Vector<Real> & u,
						Vector<Real> & u_dot,
						Vector<Real> & u_dot_dot,
						Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real delta_t_2_d2 = delta_t * delta_t / 2.;
  //  Real delta_t_d2   = delta_t / 2.;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      /// @f$ \tilde{u_{n+1}} = u_{n} +  \Delta t \dot{u}_n + (1 - 2 h \beta) \frac{\Delta t^2}{2} \ddot{u}_n @f$
      *u_val += delta_t * *u_dot_val + delta_t_2_d2 * (1 - 2 * h *  beta)* *u_dot_dot_val;

      /// @f$ \tilde{\dot{u}_{n+1}} = \dot{u}_{n} +  (1 - h \gamma) \Delta t \ddot{u}_{n} @f$
      *u_dot_val += (1 - gamma) * delta_t * *u_dot_dot_val;

      /// @f$ \tilde{\ddot{u}_{n+1}} = (1 - h ) \ddot{u}_{n} @f$
      *u_dot_dot_val *= (1 - h);
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NewmarkBeta::integrationSchemeCorrImplicit(Real delta_t,
						Vector<Real> & delta_u,
						Vector<Real> & u,
						Vector<Real> & u_dot,
						Vector<Real> & u_dot_dot,
						Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  //  Real delta_t_d2 =  delta_t / 2.;
  Real inv_beta_delta_t_2 = 1. / (beta * delta_t * delta_t);
  Real gamma_inv_beta_delta_t   = gamma / (beta * delta_t);

  Real * delta_u_val   = delta_u.values;
  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      *u_val += *delta_u_val;

      *u_dot_val += gamma_inv_beta_delta_t * *delta_u_val;

      *u_dot_dot_val += inv_beta_delta_t_2 * *delta_u_val;
    }

    delta_u_val++;
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}
