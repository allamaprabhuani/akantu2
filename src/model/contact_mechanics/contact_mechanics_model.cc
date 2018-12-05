/**
 * @file   coontact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue May 08 2012
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Contact mechanics model
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
#include "fe_engine.hh"
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */


namespace akantu {

ContactMechanicsModel::ContactMechanicsModel( Mesh & mesh, UInt dim, const ID & id,
					      const MemoryID & memory_id,
					      const ModelType model_type)
  : Model(mesh, model_type, dim, id, memory_id) {

  AKANTU_DEBUG_IN();


  //this->detector = std::make_unique<ContactDetector>(
  //	  this->mesh, id + ":contact_detector");

  AKANTU_DEBUG_OUT();
  
}


ContactMechanicsModel::~ContactMechanicsModel() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

void ContactMechanicsModel::initModel() {
  
}
  

}
