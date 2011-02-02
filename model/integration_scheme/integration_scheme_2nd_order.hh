/**
 * @file   integration_scheme_2nd_order.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug  2 12:22:05 2010
 *
 * @brief  Interface of the integrator of second order
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

  virtual ~IntegrationScheme2ndOrder() {};
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
