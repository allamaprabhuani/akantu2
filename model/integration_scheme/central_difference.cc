/**
 * @file   central_difference.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug  2 14:40:51 2010
 *
 * @brief  Implementation of the central difference scheme
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "central_difference.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

void CentralDifference::integrationSchemePred(Real delta_t,
					      Vector<Real> & u,
					      Vector<Real> & u_dot,
					      Vector<Real> & u_dot_dot,
					      Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real delta_t_2_d2 = delta_t * delta_t / 2.;
  Real delta_t_d2   = delta_t / 2.;

  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      /// @f$ u_{n+1} = u_{n} +  \Delta t * \dot{u}_n + \frac{\Delta t^2}{2} * \ddot{u}_n @f$
      *u_val += delta_t * *u_dot_val + delta_t_2_d2 * *u_dot_dot_val;

      /// @f$ \dot{u}_{n+1} = \dot{u}_{n} +  \frac{\Delta t}{2} * \ddot{u}_{n} @f$
      *u_dot_val += delta_t_d2 * *u_dot_dot_val;
    }
    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void CentralDifference::integrationSchemeCorr(Real delta_t,
					      Vector<Real> & u,
					      Vector<Real> & u_dot,
					      Vector<Real> & u_dot_dot,
					      Vector<bool> & boundary) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = u.getSize();
  UInt nb_degre_of_freedom = u.getNbComponent() * nb_nodes;

  Real delta_t_d2 =  delta_t / 2.;

  //  Real * u_val         = u.values;
  Real * u_dot_val     = u_dot.values;
  Real * u_dot_dot_val = u_dot_dot.values;
  bool * boundary_val  = boundary.values;

  for (UInt d = 0; d < nb_degre_of_freedom; d++) {
    if(!(*boundary_val)) {
      /// @f$ \dot{u}_{n+1} = \dot{u}_{n} +  \frac{\Delta t}{2} * \ddot{u}_{n+1} @f$
      *u_dot_val += delta_t_d2 * *u_dot_dot_val;
    }
    //    u_val++;
    u_dot_val++;
    u_dot_dot_val++;
    boundary_val++;
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
