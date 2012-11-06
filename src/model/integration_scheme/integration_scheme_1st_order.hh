/**
 * @file   integration_scheme_1st_order.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Jul 04 15:07:57 2011
 *
 * @brief  Interface of the time integrator of first order
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
#include "aka_common.hh"

#ifndef __AKANTU_INTEGRATION_SCHEME_1ST_ORDER_HH__
#define __AKANTU_INTEGRATION_SCHEME_1ST_ORDER_HH__

__BEGIN_AKANTU__

class IntegrationScheme1stOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  enum IntegrationSchemeCorrectorType {
    _temperature_corrector,
    _temperature_rate_corrector
  };


  virtual ~IntegrationScheme1stOrder() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:


  virtual void integrationSchemePred(Real delta_t,
				     Vector<Real> & u,
				     Vector<Real> & u_dot,
				     Vector<bool> & boundary) = 0;

  virtual void integrationSchemeCorrTemp(Real delta_t,
					 Vector<Real> & u,
					 Vector<Real> & u_dot,
					 Vector<bool> & boundary,
					 Vector<Real> & delta) = 0;

  virtual void integrationSchemeCorrTempRate(Real delta_t,
					     Vector<Real> & u,
					     Vector<Real> & u_dot,
					     Vector<bool> & boundary,
					     Vector<Real> & delta) = 0;

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

#include "generalized_trapezoidal.hh"

#endif /* __AKANTU_INTEGRATION_SCHEME_1ST_ORDER_HH__ */
