/**
 * @file   coupler_solid_contact_explicit.cc
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
#include "coupler_solid_contact.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template<SolidMechanicsModel & solid, ContactMechanicsModel & contact>
CouplerSolidContact<solid, contact>::CouplerSolidContact() {
  this->spatial_dimension =
    solid.getMesh().getSpatialDimension();
}

/* -------------------------------------------------------------------------- */
template<SolidMechancisModel & solid, ContactMechanicsModel & contact>
CouplerSolidContact<solid, contact>::~CouplerSolidContact() {

}
									     
/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & model, ContactMechanicsModel & contact>
void CouplerSolidContact<solid, contact>::solve() {
  AKANTU_DEBUG_IN();

  switch (contact.getAnalysisMethod()) {
  case _explicit_contact: {
    solid.solveStep();
    contact.solveStep();
    break;
  }
  case _implict_contact: {
    contact.solveStep();
    break;
  }
  default:
    break;
  }
  
  this->coupleModels();
  solid.solveStep();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & solid, ContactMechanicModel & contact>
void CouplerSolidContact<solid, contact>::coupleModels() {

  switch (contact.getAnalysisMethod()) {
  case _explict_contact: {
    this->coupleExternalForces();
    break;
  }
  case _implicit_contact: {
    this->coupleExternalForces();
    this->coupleStiffnessMatrices();
    break;
  }  
  default:
    break;
  }
  this->coupleExternalForces();
}
  
/* -------------------------------------------------------------------------- */
template<SolidMechancisModel & solid, ContactMechanicsModel & contact>  
void CouplerSolidContact<solid, contact>::coupleExternalForces() {

  AKANTU_TO_IMPLEMENT();

}

/* -------------------------------------------------------------------------- */
template<SolidMechancisModel & solid, ContactMechanicsModel & contact>  
void CouplerSolidContact<solid, contact>::coupleStiffnessMatrices() {

  AKANTU_TO_IMPLEMENT();
}
  
  
}
