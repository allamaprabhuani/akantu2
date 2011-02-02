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
