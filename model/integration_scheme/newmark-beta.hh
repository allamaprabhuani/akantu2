/**
 * @file   newmark-beta.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep 30 11:35:15 2010
 *
 * @brief  general case of Newmark-@f$\beta@f$
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

#ifndef __AKANTU_NEWMARK_BETA_HH__
#define __AKANTU_NEWMARK_BETA_HH__

/* -------------------------------------------------------------------------- */
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * @f$\ddot{u}_{n+1} + 2 \xi \omega \dot{u}_{n+1} + \omega^2 u_{n+1} = F_{n+1}@f$
 *
 * @f$ u_{n+1} = \tilde{u_{n+1}} + 2 \beta \frac{\Delta t^2}{2} \ddot{u}_n @f$\n
 * where @f$ \tilde{u_{n+1}} = u_{n} +  \Delta t \dot{u}_n + (1 - 2 \beta) \frac{\Delta t^2}{2} \ddot{u}_n @f$
 *
 * @f$ \dot{u}_{n+1} = \tilde{\dot{u}_{n+1}} + \gamma \frac{\Delta t}{2} * \ddot{u}_{n+1} @f$\n
 * where @f$ \tilde{\dot{u}_{n+1}} = \dot{u}_{n} +  (1 - \gamma) \Delta t \ddot{u}_{n} @f$
 *
 *
 * @f$\beta = 0@f$, @f$\gamma = \frac{1}{2}@f$ explicit central difference method\n
 * @f$\beta = \frac{1}{4}@f$, @f$\gamma = \frac{1}{2}@f$ undamped trapezoidal rule (\xi = 0)\n
 * @f$\gamma > \frac{1}{2}@f$ numerically damped integrator with damping proportional to @f$\gamma - \frac{1}{2}@f$
 *
 *
 * Stability :\n
 * Unconditionally stable for @f$\beta \geq \frac{\gamma}{2} \geq \frac{1}{4}@f$\n
 * Conditional stability:\n
 * @f$ \omega_{max} \Delta t = \frac{\xi \bar{\gamma} + \left[ \bar{\gamma} + \frac{1}{4} - \beta + \xi^2 \\bar{\gamma}^2 \right]^{\frac{1}{2}}}{\left( \frac{\gamma}{2} - \beta \right)}, \bar{\gamma} \equiv \gamma - \frac{1}{2} \geq 0 @f$
 */
class NewmarkBeta : public IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NewmarkBeta(Real beta, Real gamma) : beta(beta), gamma(gamma), h(1) {};

  ~NewmarkBeta(){};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void integrationSchemePred(Real delta_t,
			     Vector<Real> & u,
			     Vector<Real> & u_dot,
			     Vector<Real> & u_dot_dot,
			     Vector<bool> & boundary);

  void integrationSchemeCorr(Real delta_t,
			     Vector<Real> & u,
			     Vector<Real> & u_dot,
			     Vector<Real> & u_dot_dot,
			     Vector<bool> & boundary);

  void integrationSchemePredImplicit(Real delta_t,
				     Vector<Real> & u,
				     Vector<Real> & u_dot,
				     Vector<Real> & u_dot_dot,
				     Vector<bool> & boundary);

  void integrationSchemeCorrImplicit(Real delta_t,
				     Vector<Real> & delta_u,
				     Vector<Real> & u,
				     Vector<Real> & u_dot,
				     Vector<Real> & u_dot_dot,
				     Vector<bool> & boundary);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Beta, beta, Real);
  AKANTU_GET_MACRO(Gamma, gamma, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the @f$\beta@f$ parameter
  Real beta;

  /// the @f$\gamma@f$ parameter
  Real gamma;

  Real h;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "newmark-beta_inline_impl.cc"


/**
 * central difference method (explicit)
 * undamped stability condition :
 * @f$ \Delta t = \alpha \Delta t_{crit} = \frac{2}{\omega_{max}} \leq \min_{e} \frac{l_e}{c_e}
 *
 */
class CentralDifference : public NewmarkBeta {
public:
  CentralDifference() : NewmarkBeta(0, 0.5) {};
};
//#include "integration_scheme/central_difference.hh"

/// undamped trapezoidal rule (implicit)
class TrapezoidalRule : public NewmarkBeta {
public:
  TrapezoidalRule() : NewmarkBeta(0.25, 0.5) {};
};


__END_AKANTU__

#endif /* __AKANTU_NEWMARK_BETA_HH__ */
