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


/* -------------------------------------------------------------------------- */
  ContactDetection::ContactDetection(Mesh & mesh,std::string master, std::string slave, const ID & id,
				   UInt memory_id)
  : Memory(id, memory_id),
    Parsable(ParserType::_contact_detection, id),
    mesh(mesh) {
  
  this->spatial_dimension = mesh.getSpatialDimension();

  this->registerParam("master_surface", master_id, _pat_parsmod,
		      "Master surface id");
  this->registerParam("slave_surface", slave_id, _pat_parsmod,
		      "Slave surface id");

  this->master_id = master;
  this->slave_id = slave;
  
  this->getMaximalDetectionDistance();
  
}

/* -------------------------------------------------------------------------- */
void ContactDetection::getMaximalDetectionDistance() {

  AKANTU_DEBUG_IN();

  Real el_size;
  Real max_el_size = std::numeric_limits<Real>::min();

  auto & master_group =
    mesh.getElementGroup(master_id);

  for (auto & type: master_group.elementTypes(spatial_dimension - 1, _not_ghost)) {

    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    Array<Real> coord(0, nb_nodes_per_element * spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), coord, type,
					 _not_ghost);
    auto el_coord = coord.begin(spatial_dimension, nb_nodes_per_element);
    UInt nb_element = mesh.getNbElement(type);

    for (UInt el =0; el < nb_element; ++el, ++el_coord) {
      el_size = FEEngine::getElementInradius(*el_coord, type);
      max_el_size = std::max(max_el_size, el_size);
    }

    AKANTU_DEBUG_INFO("The maximum element size : "
		      << max_el_size );    
  }

  this->max_dd = max_el_size;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactDetection::search() {
  this->globalSearch();
  this->localSearch();
}
  
 
/* -------------------------------------------------------------------------- */
void ContactDetection::globalSearch() {
   
  auto & master_list =
    mesh.getElementGroup(master_id).getNodeGroup().getNodes();

  auto & slave_list =
    mesh.getElementGroup(slave_id).getNodeGroup().getNodes();
   
  BBox bbox_master(spatial_dimension);
  this->constructBoundingBox(bbox_master, master_list);
  
  BBox bbox_slave(spatial_dimension);
  this->constructBoundingBox(bbox_slave, slave_list);
  
  auto && bbox_intersection =
    bbox_master.intersection(bbox_slave);
 
  Vector<Real> center(spatial_dimension);
  bbox_intersection.getCenter(center);

  Vector<Real> spacing(spatial_dimension);
  this->computeCellSpacing(spacing);

  SpatialGrid<UInt> master_grid(spatial_dimension, spacing, center);
  this->constructGrid(master_grid, bbox_intersection, master_list);

  SpatialGrid<UInt> slave_grid(spatial_dimension, spacing, center);
  this->constructGrid(slave_grid, bbox_intersection, slave_list);
  
  if (AKANTU_DEBUG_TEST(dblDump)) {
    Mesh mesh(spatial_dimension, "save");
    master_grid.saveAsMesh(mesh);
    mesh.write("grid.msh");
  }

  AKANTU_DEBUG_INFO( "Grid Details "
		     << master_grid );
  
  // (bucket sort - construct the grid in intersected bounding box)
  // create 2 array As and Am which contains ith cell j nodes 
  // need to find nodes int he intersecting box
  // based on the position of the node, a cell can be assigned
  // given in report
  // the arrays are replaced by spatial grids 

  // TODO : increase the lower and upper bound for both bbox 
}

/* -------------------------------------------------------------------------- */
void ContactDetection::localSearch() {

  // local search
  // out of these array check each cell for closet node in that cell
  // and neighbouring cells find the actual orthogonally closet
  // check the projection of slave node on master facets connected to
  // the closet master node, if yes update the contact element with
  // slave node and master node and master surfaces connected to the
  // master node
  // these master surfaces will be needed later to update contact elements

  // PART I
  // for each cell in slave grid
  // for each node in each cell
  // check the corresponding and neighboring cell in master grid
  // find the closet master mode

  // PART II
  // Once closet master node if found ,
  // find the attached elements to it surface/line
  // mesh.getAssociatedElements( function seems to be the key)

  // PART III
  // compute the orthogonal distance from each associated element

  // PART IV
  // create contact element,
  // contact element should have one slave node and master facets

  
  
  
}

/* -------------------------------------------------------------------------- */
void ContactDetection::constructGrid(SpatialGrid<UInt> & grid, BBox & bbox,
				     const Array<UInt> & nodes_list) {

  const auto & positions = mesh.getNodes();

  auto to_grid = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s: arange(spatial_dimension)) {
      pos(s) = positions(node, s);
    }

    if (bbox.contains(pos)) {
      grid.insert(node, pos);
    }
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_grid);
}

/* -------------------------------------------------------------------------- */  
void ContactDetection::constructBoundingBox(BBox & bbox, const Array<UInt> & nodes_list) {

  const auto & positions = mesh.getNodes();
  
  auto to_bbox = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s: arange(spatial_dimension)) {
      pos(s) = positions(node, s);
    }
    
    bbox += pos;
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_bbox); 
}

/* -------------------------------------------------------------------------- */
void ContactDetection::computeCellSpacing(Vector<Real> & spacing) {

  for (UInt s: arange(spatial_dimension)) 
    spacing(s) = std::sqrt(2.0) * max_dd;
    
}
  

  
} // akantu
