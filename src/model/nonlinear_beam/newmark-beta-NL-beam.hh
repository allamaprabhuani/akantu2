/**
 * @file   newmark-beta.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Sat Sep 12 2020
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "newmark-beta.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NEWMARK_BETA_NL_BEEAM_HH_
#define AKANTU_NEWMARK_BETA_NL_BEAM_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {


class NewmarkBetaNLBeam : public CentralDifference {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NewmarkBetaNLBeam(DOFManager & dof_manager, const ID & dof_id);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                 Array<Real> & u_dot_dot,
                 const Array<bool> & blocked_dofs) const override;

  void corrector(Real delta_t, Array<Real> & u,
                 Array<Real> & u_dot, Array<Real> & u_dot_dot,
                 const Array<bool> & blocked_dofs) const;


 
};

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_NEWMARK_BETA_NL_BEAM_HH_ */
