/**
 * @file   central_difference.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug  2 14:22:56 2010
 *
 * @brief explicit hyperbolic 2nd order central difference
 *        Newmark-beta with gama = 1/2 and beta = 0
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CENTRAL_DIFFERENCE_HH__
#define __AKANTU_CENTRAL_DIFFERENCE_HH__

__BEGIN_AKANTU__

class CentralDifference : public IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

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

#endif /* __AKANTU_CENTRAL_DIFFERENCE_HH__ */
