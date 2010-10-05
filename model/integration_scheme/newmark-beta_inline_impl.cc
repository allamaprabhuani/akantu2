/**
 * @file   newmark-beta_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep 30 11:45:20 2010
 *
 * @brief  implementation of the newmark-@f$\beta@f$ integration scheme
 *
 * @section LICENSE
 *
 * <insert license here>
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
    //    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}
