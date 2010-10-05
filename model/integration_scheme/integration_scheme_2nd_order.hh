/**
 * @file   integration_scheme_2nd_order.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug  2 12:22:05 2010
 *
 * @brief  Interface of the integrator of second order
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__
#define __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void integrationSchemePred(Real delta_t,
				     Vector<Real> & u,
				     Vector<Real> & u_dot,
				     Vector<Real> & u_dot_dot,
				     Vector<bool> & boundary) = 0;

  virtual void integrationSchemeCorr(Real delta_t,
				     Vector<Real> & u,
				     Vector<Real> & u_dot,
				     Vector<Real> & u_dot_dot,
				     Vector<bool> & boundary) = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};

__END_AKANTU__

#include "integration_scheme/newmark-beta.hh"

#endif /* __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__ */
