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

  
  ContactDetection::ContactDetection(Mesh & mesh)
    : mesh(mesh) {
    
    this->spatial_dimension = mesh.getSpatialDimension();
}

  
/* -------------------------------------------------------------------------- */
void ContactDetection::globalSearch() {
  // global search
  // bounding box to determine possible slave nodes and master nodes
  // create bounding boxes from slave and master surfaces
  
  auto & master_list =
    mesh.getElementgroup(master_id).getNodeGroup().getNodes();

  auto & slave_list =
    mesh.getElementgroup(slave_id).getNodeGroup().getNodes();
   
  BBox bbox_master(spatial_dimension);
  this->constructBoundingBox(bbox_master, master_list);
  
  BBox bbox_slave(spatial_dimension);
  this->constructBoundingBox(bbox_slave, slave_list);
  
  auto && bbox_intersection =
    bbox_master.intersection(bbox_slave);
  
  // (bucket sort - construct the grid in intersected bounding box)
  // create 2 array As and Am which contains ith cell j nodes 
  // need to find nodes int he intersecting box
  // based on the position of the node, a cell can be assigned
  // given in report

   
  // get the upper bound for bbox
  // get the lower bound for bbox
  // compute the center of the box
  const auto & lower = bbox_intersection.getLowerBounds();
  const auto & upper = bbox_intersection.getUpperBounds();

  Vector<Real> center = upper;
  center += lower;
  center /= 2.;

  // define the spacing or size of cells
  //  SpatialGrid<Element> grid(spatial_dimension, spacing, center);
  // 
  
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

/* -------------------------------------------------------------------------- */
void ContactDetection::constructGrid(SpatialGrid<Element> & grid) {
  Vector<Real> bary(spatial_dimension);
  Element el;
  el.ghost_type = _not_ghost;

  auto it = mesh.firstType(spatial_dimension);
  auto last_type = mesh.lastType(spatial_dimension);
  for (; it != last_type; ++it) {
    UInt nb_element = mesh.getNbElement(*it);
    el.type = *it;
    for (UInt e = 0; e < nb_element; ++e) {
      el.element = e;
      mesh.getBaryCenter(el, bary);
      grid.insert(el, bary);
    }
  }
}

/* -------------------------------------------------------------------------- */  
  void ContactDetection::constructBoundingBox(BBox & bbox, const Array<UInt> & nodes_list) {

  auto to_position = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (Uint s: arange(spatial_dimension)) {
      pos(s) = positions(node, s);
    }
    auto && info = NodeInfo(node, pos, direction);
    bbox += info.position;
    return info;
  };

  std::vector<NodeInfo> nodes(master_list.size());

  std::transform(nodes_list.begin(), nodes_list.end(), nodes.begin(),
		 to_position);
  
}

  
} // akantu
