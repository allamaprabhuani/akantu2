/**
 * @file  contact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Sep 10 2018
 * @date last modification: Mon Sep 10 2018
 *
 * @brief  Model of Contact Mechanics
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
#include "contact_mechanics_model.hh"

namespace akantu {

  ContactMechanicsModel::ContactMechanicsModel(Mesh & mesh, UInt dim, const ID & id,
					       const MemoryID & memory_id,
					       const ModelType model_type) 
    : Model(mesh, model_type, dim, id,  memory_id),
      mesh(mesh) {

  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::initFullImpl(const ModelOptions & options) {
  AKANTU_DEBUG_IN();

  /*resolution_method = options.analysis_method;
  if (!this->hasDefaultResolution()) {
    this->initNewResolution(this->resolution_method);
    }*/
  
  
  AKANTU_DEBUG_OUT();

}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::createContactElements() {

  // based on type of detection object
  // pass the surfaces master and slave to create contact pairs
  // let say detection system is node to segement and bucket sort

  // only consider those slave nodes which doesnot have any othogonal
  // projection on its contact element
  
  // global search
  // bounding box to determine possible slave nodes and master nodes
  // (bucket sort - construct the grid in intersected bounding box)
  // create 2 array As and Am which contains ith cell j nodes 
  // 

  // local search
  // out of these arraye check each cell for closet node in that cell
  // and neighbouring cells find the actual orthogonally closet
  // check the projection of slave node on master facets connected to
  // the closet master node, if yes update teh contact element with
  // slave node and master node and master surfaces connected to the
  // master node
  // these master surfaces will be needed later to update contact elements
  
  // once find insert it in contact_elements
  
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsModel::updateContactElements() {


}

/* -------------------------------------------------------------------------- */
//void ContactMechanicsModel::initResolution() {

//}

/* -------------------------------------------------------------------------- */
//void ContactMechanicsModel::solve() {

//}

/* -------------------------------------------------------------------------- */  
void ContactMechanicsModel::printself(std::ostream & stream, int indent) const {
  /*std::string space;
  for (Int i = 0; i < indent; space += AKANTU_INDENT) 
    ;

  stream << space << "Contact Mechanics Model [" << std::endl;
  stream << space << " + id                  : " << id << std::endl;
  stream << space << " + spatial dimension   : " << model.getSpatialDimension()
	 << std::endl;
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;*/
}
  
} // namespace akantu
