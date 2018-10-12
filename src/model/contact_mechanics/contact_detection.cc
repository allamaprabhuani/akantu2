/**
 * @file contact_detection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Sep 12 2018
 * @date last modification: Fri Sep 21 2018
 *
 * @brief  Mother class for all detection algorithms
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
#include "contact_detection.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

ContactDetection::ContactDetection() {
  

  
}

/* -------------------------------------------------------------------------- */
void ContactDetection::globalSearch() {
  // global search
  // bounding box to determine possible slave nodes and master nodes
  // create bounding boxes from slave and master surfaces
  // BBbox 

  // (bucket sort - construct the grid in intersected bounding box)
  // create 2 array As and Am which contains ith cell j nodes 
  // need to find nodes int he intersecting box
  // based on the position of the node, a cell can be assigned
  // given in report
  
}

/* -------------------------------------------------------------------------- */
void ContactDetection::localSearch() {

  // local search
  // out of these arraye check each cell for closet node in that cell
  // and neighbouring cells find the actual orthogonally closet
  // check the projection of slave node on master facets connected to
  // the closet master node, if yes update teh contact element with
  // slave node and master node and master surfaces connected to the
  // master node
  // these master surfaces will be needed later to update contact elements

}
  


} // akantu
