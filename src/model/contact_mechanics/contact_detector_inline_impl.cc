/**
 * @file contact_detection.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Apr 29 2019
 * @date last modification: Mon Apr 29 2019
 *
 * @brief  inine implementation of the contact detector class
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
#include "contact_detector.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__
#define __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool ContactDetector::checkValidityOfProjection(Vector<Real> & projection) {

  UInt nb_xi_inside = 0;
  Real epsilon = 1e-3;

  for (auto xi : projection) {
    if (xi >= -1.0 - epsilon and xi <= 1.0 + epsilon) 
      nb_xi_inside++;
  }  

  if (nb_xi_inside == projection.size()) 
    return true;
  
  return false;
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::coordinatesOfElement(const Element & el, Matrix<Real> & coords) {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = mesh.getConnectivity(el.type, _not_ghost)
    .begin(nb_nodes_per_element)[el.element]; 

  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s: arange(spatial_dimension)) {
      coords(s, n) = this->positions(node, s);
    }
  }
}
  
/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeCellSpacing(Vector<Real> & spacing) {
  for (UInt s: arange(spatial_dimension)) 
    spacing(s) = std::sqrt(2.0) * max_dd;
}

/* -------------------------------------------------------------------------- */  
inline void ContactDetector::constructBoundingBox(BBox & bbox, const Array<UInt> & nodes_list) {
  
  auto to_bbox = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s: arange(spatial_dimension)) {
      pos(s)  = this->positions(node, s);
    }
    bbox += pos;
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_bbox);

  auto & lower_bound = bbox.getLowerBounds();
  auto & upper_bound = bbox.getUpperBounds();

  for (UInt s: arange(spatial_dimension)) {
    lower_bound(s) -= this->max_bb;
    upper_bound(s) += this->max_bb;
  }

  std::cerr << bbox << std::endl;
  
  AKANTU_DEBUG_INFO("BBox" << bbox);
}


/* -------------------------------------------------------------------------- */
inline void ContactDetector::constructGrid(SpatialGrid<UInt> & grid, BBox & bbox,
					   const Array<UInt> & nodes_list) {
  auto to_grid = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s: arange(spatial_dimension)) {
      pos(s) = this->positions(node, s);
    }

    if (bbox.contains(pos)) {
      grid.insert(node, pos);
    }
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_grid);
}
  
/* -------------------------------------------------------------------------- */
template<Surface id>
inline std::string ContactDetector::getSurfaceId() {
  return surfaces[id];
}

/* -------------------------------------------------------------------------- */
template<Surface id>
inline void ContactDetector::setSurfaceId(const std::string name) {
  surfaces[id] = name;
}
  
/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeMaximalDetectionDistance() {

  AKANTU_DEBUG_IN();

  Real elem_size;
  Real max_elem_size = std::numeric_limits<Real>::min();

  auto & master_group =
    mesh.getElementGroup(surfaces[Surface::master]);

  for (auto type:
	 master_group.elementTypes(spatial_dimension - 1, _not_ghost, _ek_regular)) {
    
    const auto & element_ids = master_group.getElements(type);    
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
   
    Element elem;
    elem.type = type;
    for (auto el : element_ids) {
      elem.element = el;
      Matrix<Real> elem_coords(spatial_dimension, nb_nodes_per_element);
      this->coordinatesOfElement(elem, elem_coords);

      elem_size = FEEngine::getElementInradius(elem_coords, type);
      max_elem_size = std::max(max_elem_size, elem_size);
    }

    AKANTU_DEBUG_INFO("The maximum element size : "
		      << max_elem_size );
  }

  this->max_dd = max_elem_size;
  this->max_bb = max_elem_size;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt> ContactDetector::constructConnectivity(UInt & slave, const Element & master) {
    
  Vector<UInt> master_conn = this->mesh.getConnectivity(master);

  Vector<UInt> elem_conn(master_conn.size() + 1);

  elem_conn[0] = slave;
  for (UInt i = 1; i < elem_conn.size(); ++i) {
    elem_conn[i] = master_conn[i-1];
  }    

  return elem_conn;
}
  
/* -------------------------------------------------------------------------- */  

} //akantu

#endif /*  __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__ */
