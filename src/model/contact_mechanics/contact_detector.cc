/**
 * @file contact_detector.cc
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
#include "contact_detector.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
 ContactDetector::ContactDetector(Mesh & mesh, const ID & id,  UInt memory_id)
   : ContactDetector(mesh, mesh.getNodes(), id, memory_id) {
  
 }
  
/* -------------------------------------------------------------------------- */
ContactDetector::ContactDetector(Mesh & mesh, Array<Real> & positions, const ID & id, UInt memory_id)
  : Memory(id, memory_id),
    Parsable(ParserType::_contact_detector, id),
    mesh(mesh),
    positions(positions) {

  AKANTU_DEBUG_IN();  

  this->spatial_dimension = mesh.getSpatialDimension();
    
  this->mesh.fillNodesToElements(this->spatial_dimension - 1);  

  this->parseSection();

  this->computeMaximalDetectionDistance();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactDetector::parseSection() {

  const Parser & parser = getStaticParser();

  const ParserSection & section =
    *(parser.getSubSections(ParserType::_contact_detector).first);

  auto type = section.getParameterValue<std::string>("type");
  if (type == "implicit") {
    this->detection_type = _implicit;
  }
  else if (type == "explicit"){
    this->detection_type = _explicit;
  }
  else {
    AKANTU_ERROR("Unknown detection type : " << type);
  }
  
  surfaces[Surface::master] = section.getParameterValue<std::string>("master");
  surfaces[Surface::slave ] = section.getParameterValue<std::string>("slave");
}
     
/* -------------------------------------------------------------------------- */
void ContactDetector::search(std::map<UInt, ContactElement> & contact_map) {
  
  SpatialGrid<UInt> master_grid(spatial_dimension);

  SpatialGrid<UInt> slave_grid(spatial_dimension);
 
  this->globalSearch(slave_grid, master_grid);

  this->localSearch(slave_grid, master_grid);

  this->constructContactMap(contact_map);
}
   
/* -------------------------------------------------------------------------- */
void ContactDetector::globalSearch(SpatialGrid<UInt> & slave_grid,
				   SpatialGrid<UInt> & master_grid) {
  
  auto & master_list =
    mesh.getElementGroup(surfaces[Surface::master]).getNodeGroup().getNodes();

  auto & slave_list =
    mesh.getElementGroup(surfaces[Surface::slave]).getNodeGroup().getNodes();
   
  BBox bbox_master(spatial_dimension);
  this->constructBoundingBox(bbox_master, master_list);
  
  BBox bbox_slave(spatial_dimension);
  this->constructBoundingBox(bbox_slave, slave_list);
  
  auto && bbox_intersection =
    bbox_master.intersection(bbox_slave);
  
  AKANTU_DEBUG_INFO( "Intersection BBox "
		     << bbox_intersection );
 
  Vector<Real> center(spatial_dimension);
  bbox_intersection.getCenter(center);

  Vector<Real> spacing(spatial_dimension);
  this->computeCellSpacing(spacing);
 
  master_grid.setCenter(center);
  master_grid.setSpacing(spacing);
  this->constructGrid(master_grid, bbox_intersection, master_list);

  slave_grid.setCenter(center);
  slave_grid.setSpacing(spacing);
  this->constructGrid(slave_grid, bbox_intersection, slave_list);
    
  // search slave grid nodes in contactelement array and if they exits
  // and still have orthogonal projection on its associated master
  // facetremove it from the spatial grid or do not consider it for
  // local search, maybe better option will be to have spatial grid of
  // type node info and one of the variable of node info should be
  // facet already exits
  // so contact elements will be updated based on the above
  // consideration , this means only those contact elements will be
  // keep whose slave node is still in intersection bbox and still has
  // projection in its master facet
  // also if slave node is already exists in contact element and
  // orthogonal projection does not exits then search the associated
  // master facets with the current master facets within a given
  // radius , this is subjected to computational cost as searching
  // neighbbor cells can be more effective.
}

/* -------------------------------------------------------------------------- */
void ContactDetector::localSearch(SpatialGrid<UInt> & slave_grid,
				  SpatialGrid<UInt> & master_grid) {

  // local search
  // out of these array check each cell for closet node in that cell
  // and neighbouring cells find the actual orthogonally closet
  // check the projection of slave node on master facets connected to
  // the closet master node, if yes update the contact element with
  // slave node and master node and master surfaces connected to the
  // master node
  // these master surfaces will be needed later to update contact
  // elements

  //Array<UInt> slave_nodes;
  //Array<UInt> master_nodes;

  /*BBox bbox_master_grid(spatial_dimension); 
  BBox bbox_slave_grid(spatial_dimension);

  auto create_bbox = [&](auto & grid, auto & bbox) {
    auto upper_bound = grid.getUpperBounds();
    auto lower_bound = grid.getLowerBounds();
    for (UInt s: arange(spatial_dimension)) {
      lower_bound(s) -= this->max_bb;
      upper_bound(s) += this->max_bb;
    }

    bbox += lower_bound;
    bbox += upper_bound;
  };

  create_bbox(master_grid, bbox_master_grid);
  create_bbox(slave_grid, bbox_slave_grid);
  
  auto && bbox_intersection =
  bbox_master_grid.intersection(bbox_slave_grid);*/

  contact_pairs.clear();
  
  /// find the closet master node for each slave node
  for (auto && cell_id : slave_grid) {
    /// loop over all the slave nodes of the current cell    
    for (auto && slave_node: slave_grid.getCell(cell_id)) {

      bool pair_exists = false;
      
      Vector<Real> pos(spatial_dimension);
      for (UInt s: arange(spatial_dimension)) 
	pos(s) = this->positions(slave_node, s);
            
      //if (!bbox_intersection.contains(pos)) {
      //	continue;
      //}

      Real closet_distance = std::numeric_limits<Real>::max();
      UInt closet_master_node;

      /// loop over all the neighboring cells of the current cell
      for (auto && neighbor_cell : cell_id.neighbors()) {
	/// loop over the data of neighboring cells from master grid	
	for (auto && master_node : master_grid.getCell(neighbor_cell)) {
	  
	  Vector<Real> pos2(spatial_dimension);
	  for (UInt s: arange(spatial_dimension)) 
	    pos2(s) = this->positions(master_node, s);
	  	  
	  Real distance = pos.distance(pos2);
  
	  if (distance <= closet_distance) {
	    closet_master_node = master_node;
	    closet_distance = distance;
	    pair_exists = true;
	  }
	}	
      }

      if (pair_exists) {
	contact_pairs.push_back(slave_node);
	contact_pairs.push_back(closet_master_node);
      }
      
    }
  }  
}

/* -------------------------------------------------------------------------- */
void ContactDetector::constructContactMap(std::map<UInt, ContactElement> & contact_map) {

  auto get_index = [&](auto & gaps, auto & projections) {

    UInt index;
    Real gap_min  = std::numeric_limits<Real>::max();

    UInt counter = 0;
    for (auto && values : zip(gaps,
			      make_view(projections, spatial_dimension -1))) {
      auto & gap        = std::get<0>(values);
      auto & projection = std::get<1>(values);
      
      bool is_valid = this->checkValidityOfProjection(projection);
      
      if (is_valid and gap <= gap_min) {
	gap_min = gap;
	index = counter;
      }
      counter++;
    }

    /// TODO: adhoc fix to assign a master element in case the
    /// projection does not lie in the extended element. As it is
    /// tolerance based
    if (index >= gaps.size()) {
      auto gap_min_it = std::min_element(gaps.begin(), gaps.end());
      auto index_it = std::find(gaps.begin(), gaps.end(), *gap_min_it);
      index  = *index_it;
    }
    
    return index;
  };


  auto get_connectivity = [&](auto & slave, auto & master) {
    Vector<UInt> master_conn = this->mesh.getConnectivity(master);

    Vector<UInt> elem_conn(master_conn.size() + 1);

    elem_conn[0] = slave;
    for (UInt i = 1; i < elem_conn.size(); ++i) {
      elem_conn[i] = master_conn[i-1];
    }    

    return elem_conn;
  };

  
  for (auto && pairs : make_view(contact_pairs, 2)) {

    const auto & slave_node  = pairs(0);
    const auto & master_node = pairs(1);

    Array<Element> elements;
    this->mesh.getAssociatedElements(master_node, elements);

    Array<Real>    gaps(elements.size(),        1,                    "gaps");
    Array<Real>    normals(elements.size(),     spatial_dimension,    "normals");
    Array<Real>    projections(elements.size(), spatial_dimension -1, "projections");

    this->computeOrthogonalProjection(slave_node, elements,
				      normals, gaps, projections);
    
    auto index = get_index(gaps, projections);

    auto connectivity = get_connectivity(slave_node, elements[index]);
    
    contact_map[slave_node].setMaster(elements[index]);
    contact_map[slave_node].setGap(gaps[index]);
    contact_map[slave_node].setNormal(Vector<Real>(normals.begin(spatial_dimension)[index], true));
    contact_map[slave_node].setProjection(Vector<Real>(projections.begin(spatial_dimension - 1)[index], true));
    contact_map[slave_node].setConnectivity(connectivity);

    /// number of surface tangents will be equal to dimension of
    /// surface i.e. spatial_dimension - 1 and nb of components will
    /// still be equal to spatial dimension  
    Matrix<Real> tangents(spatial_dimension - 1, spatial_dimension);
    this->computeTangentsOnElement(contact_map[slave_node].master,
				   contact_map[slave_node].projection,
				   tangents);
    contact_map[slave_node].setTangent(tangents);
  }
  
}
  
/* -------------------------------------------------------------------------- */
void ContactDetector::computeOrthogonalProjection(const UInt & node,
						  const Array<Element> & elements,
						  Array<Real> & normals, Array<Real> & gaps,
						  Array<Real> & projections) {

  Vector<Real> query(spatial_dimension);
  for (UInt s: arange(spatial_dimension)) 
    query(s) = this->positions(node, s);
  
  for (auto && values :
	 zip( elements,
	      gaps, 
	      make_view(normals , spatial_dimension),
	      make_view(projections, spatial_dimension - 1))) {
    
    const auto & element = std::get<0>(values);
    auto & gap           = std::get<1>(values);
    auto & normal        = std::get<2>(values);
    auto & projection    = std::get<3>(values);
    
    this->computeNormalOnElement(element, normal);
   
    Vector<Real> real_projection(spatial_dimension);
    this->computeProjectionOnElement(element, normal, query,
				     projection, real_projection);     
    
    Vector<Real> distance(spatial_dimension);
    distance = query - real_projection;
    gap = Math::norm(spatial_dimension, distance.storage());
    
    Vector<Real> direction = distance.normalize(); 
    Real cos_angle = direction.dot(normal);

    Real tolerance = 1e-8;

    /// TODO: adhoc fix to ensure that normal is always into the slave
    /// surface. However, it doesnot work if gap is 0 as cos angle is
    /// a nan value
    if (std::abs(cos_angle + 1) <= tolerance) {
      normal *= -1.0;
    }
    
    if (std::abs(cos_angle - 1) <= tolerance and detection_type == _explicit) {
      gap *= -1;
    }
  }

}
 
/* -------------------------------------------------------------------------- */
void ContactDetector::computeNormalOnElement(const Element & element, Vector<Real> & normal) {
  
  Matrix<Real> vectors(spatial_dimension, spatial_dimension - 1);
  this->vectorsAlongElement(element, vectors);
 
  switch (this->spatial_dimension) {
  case 2: {
    Math::normal2(vectors.storage(), normal.storage());
    break;
  }
  case 3: {
    Math::normal3(vectors(0).storage(), vectors(1).storage(), normal.storage());
    break;
  }  
  default: { AKANTU_ERROR("Unknown dimension : " << spatial_dimension); }
  }
        
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeProjectionOnElement(const Element & element,
						 const Vector<Real> & normal,
						 const Vector<Real> & query,
						 Vector<Real> & natural_projection,
						 Vector<Real> & real_projection) {
  
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  
  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  this->coordinatesOfElement(element, coords);

  Vector<Real> point(coords(0));  
  Real alpha = (query - point).dot(normal);
  real_projection = query - alpha * normal;
  
  this->computeNaturalProjection(element, real_projection, natural_projection);  
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeNaturalProjection(const Element & element, Vector<Real> & real_projection,
					       Vector<Real> & natural_projection) {

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type,
					 _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(this->positions, nodes_coord.storage(),
				     elem_val + element.element * nb_nodes_per_element,
				     nb_nodes_per_element, spatial_dimension);
  
#define GET_NATURAL_COORDINATE(type)					\
  ElementClass<type>::inverseMap(real_projection, nodes_coord, natural_projection) 
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NATURAL_COORDINATE);
#undef GET_NATURAL_COORDINATE
}

  
/* -------------------------------------------------------------------------- */
void ContactDetector::vectorsAlongElement(const Element & el, Matrix<Real> & vectors) {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  this->coordinatesOfElement(el, coords);

  switch (spatial_dimension) {
  case 2: {
    vectors(0) = Vector<Real>(coords(1)) - Vector<Real>(coords(0));
    break;
  }
  case 3: {
    vectors(0) = Vector<Real>(coords(1)) - Vector<Real>(coords(0));
    vectors(1) = Vector<Real>(coords(2)) - Vector<Real>(coords(0));
    break;
  } 
  default: { AKANTU_ERROR("Unknown dimension : " << spatial_dimension); }
  }
  
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeTangentsOnElement(const Element & el, Vector<Real> & projection, Matrix<Real> & tangents) {

  const ElementType & type  = el.type;
  
  UInt nb_nodes_master = Mesh::getNbNodesPerElement(type);

  Vector<Real> shapes(nb_nodes_master);
  Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);
       
#define GET_SHAPES_NATURAL(type)				\
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL  

#define GET_SHAPE_DERIVATIVES_NATURAL(type)				\
  ElementClass<type>::computeDNDS(projection, shapes_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL
    
  
  Matrix<Real> coords(spatial_dimension, nb_nodes_master);
  coordinatesOfElement(el, coords);
  
  tangents.mul<false, true>(shapes_derivatives, coords);

  auto temp_tangents = tangents.transpose();
  for (UInt i = 0; i < spatial_dimension -1; ++i) {
    auto temp = Vector<Real>(temp_tangents(i));
    temp_tangents(i) = temp.normalize();
  }

  tangents = temp_tangents.transpose();
}
 
/* -------------------------------------------------------------------------- */
void ContactDetector::normalProjection(const Element & /*el*/, const Vector<Real> & /*slave_coord*/,
				       Vector<Real> & /*natural_coord*/, Real & /*tolerance*/) {

  /*Real fmin;

  
  auto update_fmin = [&fmin, &slave_coord, &node_coords, &natural_coord]() {
    Vector<Real> physical_guess_v(physical_guess.storage(), spatial_dimension);
    // interpolate on natual coordinate and get physical guess
    // compute gradient or jacobian on natural cooordiante
    // f = slave_coord - physical_guess;
    
    return fmin;
  };

  auto closet_point_error = update_fmin();

  while (tolerance < closet_point_error) {

    // compute gradiend on natural coordinate
    // compute second variation of shape function at natural coord
    
    
    closet_point_error = update_fmin();
    }*/
  
}
  
   
  
} // akantu
