/**
 * @file   coupler_solid_contact_explicit.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu Jan 17 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model in explicit 
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_commmon.hh"
#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COUPLER_SOLID_CONTACT_HH__
#define __AKANTU_COUPLER_SOLID_CONTACT_HH__

namespace akantu {

template<SolidMechanicsModel solid, ContactMechanicsModel contact>
class CouplerSolidContact  {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructor                                                   */
  /* ------------------------------------------------------------------------ */
public:
  CouplerSolidContact();

  ~CouplerSolidContact();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the coupled model
  void solve();

private:
  /// does the appropriate coupling between models
  void coupleModels();

  /// couples external forces between models
  void coupleExternalForces();

  /// couples stiffness matrices between models
  void coupleStiffnessMatrices();

};

} // akantu

#endif /* __COUPLER_SOLID_CONTACT_HH__  */
