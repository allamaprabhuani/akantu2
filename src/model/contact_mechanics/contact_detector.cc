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
  this->getMaximalDetectionDistance();
  
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
  
  surfaces[Surface::master] = section.getParameterValue<std::string>("master_surface");
  surfaces[Surface::slave ] = section.getParameterValue<std::string>("slave_surface");
}
    
/* -------------------------------------------------------------------------- */
void ContactDetector::getMaximalDetectionDistance() {

  AKANTU_DEBUG_IN();

  Real el_size;
  Real max_el_size = std::numeric_limits<Real>::min();

  auto & master_group =
    mesh.getElementGroup(surfaces[Surface::master]);

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
  this->max_bb = max_el_size;
  
  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
  void ContactDetector::search(std::map<UInt, ContactElement> & contact_map) {
  this->globalSearch(contact_map);
}
  
 
/* -------------------------------------------------------------------------- */
void ContactDetector::globalSearch(std::map<UInt, ContactElement> & contact_map) {
  
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
    
  auto & master_surface_list =
    mesh.getElementGroup(surfaces[Surface::master]).getNodeGroup().getNodes();

  auto & slave_surface_list =
    mesh.getElementGroup(surfaces[Surface::slave]).getNodeGroup().getNodes();
 
  SpatialGrid<UInt> master_grid(spatial_dimension, spacing, center);
  this->constructGrid(master_grid, bbox_intersection, master_surface_list);

  SpatialGrid<UInt> slave_grid(spatial_dimension, spacing, center);
  this->constructGrid(slave_grid, bbox_intersection, slave_surface_list);
  
  if (AKANTU_DEBUG_TEST(dblDump)) {
    Mesh mesh(spatial_dimension, "save");
    master_grid.saveAsMesh(mesh);
    mesh.write("master_grid.msh");
  }

  if (AKANTU_DEBUG_TEST(dblDump)) {
    Mesh mesh2(spatial_dimension, "save");
    slave_grid.saveAsMesh(mesh2);
    mesh2.write("slave_grid.msh");
  }
  
  AKANTU_DEBUG_INFO( "Grid Details " << master_grid );
  // search slave grid nodes in contactelement array and if they exits
  // and still have orthogonal projection on its associated master
  // facetremove it from the spatial grid or do not consider it for
  // local search, maybe better option will be to have spatial grid of
  // type node info and one of the variable of node info should be
  // facet already exits
  // so contact eleemnts will be updated based on the above
  // consideration , this means only those contact elements will be
  // keep whose slave node is still in intersection bbox and still has
  // projection in its master facet
  // also if slave node is already exists in contact element and
  // orthogonal projection does not exits then search the associated
  // master facets with the current master facets within a given
  // radius , this is subjected to computational cost as searching
  // neighbbor cells can be more effective.
  this->localSearch(slave_grid, master_grid, contact_map);
}

/* -------------------------------------------------------------------------- */
void ContactDetector::localSearch(SpatialGrid<UInt> & slave_grid,
				  SpatialGrid<UInt> & master_grid,
				  std::map<UInt, ContactElement> & contact_map) {

  // local search
  // out of these array check each cell for closet node in that cell
  // and neighbouring cells find the actual orthogonally closet
  // check the projection of slave node on master facets connected to
  // the closet master node, if yes update the contact element with
  // slave node and master node and master surfaces connected to the
  // master node
  // these master surfaces will be needed later to update contact
  // elements

  Array<UInt> slave_nodes;
  Array<UInt> master_nodes;

  BBox bbox_master_grid(spatial_dimension); 
  BBox bbox_slave_grid(spatial_dimension);

  bbox_master_grid += master_grid.getUpperBounds();
  bbox_master_grid += master_grid.getLowerBounds();

  bbox_slave_grid += slave_grid.getUpperBounds();
  bbox_slave_grid += slave_grid.getLowerBounds();

  auto && bbox_intersection =
    bbox_master_grid.intersection(bbox_slave_grid);
  
  // find the closet master node for each slave node
  for (auto && cell_id : slave_grid) {
    AKANTU_DEBUG_INFO("Looping on next cell");
    
    for (auto && q1: slave_grid.getCell(cell_id)) {
       
      Vector<Real> pos(spatial_dimension);
      for (UInt s: arange(spatial_dimension)) {
	pos(s) = this->positions(q1, s);
      }

      if (!bbox_intersection.contains(pos)) {
	continue;
      }

      Real closet_distance = std::numeric_limits<Real>::max();
      UInt closet_master_node;
     
      // loop over all the neighboring cells of the current cells
      for (auto && neighbor_cell : cell_id.neighbors()) {

	// loop over the data of neighboring cells from master grid	
	for (auto && q2 : master_grid.getCell(neighbor_cell)) {
	  
	  AKANTU_DEBUG_INFO("Looping on neighbor cell in master");
	  Vector<Real> pos2(spatial_dimension);
	  for (UInt s: arange(spatial_dimension)) {
	    pos2(s) = this->positions(q2, s);
	  }

	  Real distance = pos.distance(pos2);
  
	  if (distance <= closet_distance) {
	    closet_master_node = q2;
	    closet_distance = distance;
	  }
	}	
      }
	
      slave_nodes.push_back(q1);
      master_nodes.push_back(closet_master_node);
    }
  }  
  
  for (auto && values : zip(slave_nodes, master_nodes)) {
    const auto & slave_node  = std::get<0>(values);
    const auto & master_node = std::get<1>(values);

    Array<Element> elements;
    this->mesh.getAssociatedElements(master_node, elements);
    
    auto normals  = std::make_unique<Array<Real>>(elements.size(),
						  spatial_dimension, "normals");
    auto gaps     = std::make_unique<Array<Real>>(elements.size(),
						  1,                 "gaps"); 
    auto natural_projections = std::make_unique<Array<Real>>(elements.size(),
							     spatial_dimension - 1, "projections");
    
    this->computeOrthogonalProjection(slave_node, elements,
				      *normals, *gaps, *natural_projections);
      
    auto minimum = std::min_element( gaps->begin(), gaps->end());
    auto index   = std::distance(    gaps->begin(), minimum);
       
    Vector<UInt> master_conn =
      this->mesh.getConnectivity(elements[index]);

    Vector<UInt> elem_conn(master_conn.size() + 1);
    elem_conn[0] = slave_node;
    for (UInt i = 1; i < elem_conn.size(); ++i) {
      elem_conn[i] = master_conn[i-1];
    }    
   
    contact_map[slave_node] = ContactElement(elements[index]);
    contact_map[slave_node].gap = (*gaps)[index];
    contact_map[slave_node].normal =
      Vector<Real>(normals->begin(spatial_dimension)[index], true);
    contact_map[slave_node].projection =
      Vector<Real>(natural_projections->begin(spatial_dimension - 1)[index], true);
    contact_map[slave_node].connectivity = elem_conn;    
  }
}

  
/* -------------------------------------------------------------------------- */
void ContactDetector::constructGrid(SpatialGrid<UInt> & grid, BBox & bbox,
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
void ContactDetector::constructBoundingBox(BBox & bbox, const Array<UInt> & nodes_list) {
  
  auto to_bbox = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s: arange(spatial_dimension)) {
      pos(s)  = this->positions(node, s);
    }
    
    bbox += pos;
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_bbox);

  AKANTU_DEBUG_INFO("BBox" << bbox);
  
  auto & lower_bound = bbox.getLowerBounds();
  auto & upper_bound = bbox.getUpperBounds();

  for (UInt s: arange(spatial_dimension)) {
    lower_bound(s) -= this->max_bb;
    upper_bound(s) += this->max_bb;
  }
  
  AKANTU_DEBUG_INFO("BBox" << bbox);
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeCellSpacing(Vector<Real> & spacing) {

  for (UInt s: arange(spatial_dimension)) 
    spacing(s) = std::sqrt(2.0) * max_dd;
    
}
  
/* -------------------------------------------------------------------------- */
void ContactDetector::computeOrthogonalProjection(const UInt & node,
						  const Array<Element> & elements,
						  Array<Real> & normals, Array<Real> & gaps,
						  Array<Real> & natural_projections) {

  Vector<Real> query(spatial_dimension);
  for (UInt s: arange(spatial_dimension)) {
    query(s) = this->positions(node, s);
  }

  for (auto && values :
	 zip( elements,
	      gaps,
	      make_view(normals , spatial_dimension),
	      make_view(natural_projections, spatial_dimension - 1))) {
    const auto & element = std::get<0>(values);
    auto & gap           = std::get<1>(values);
    auto & normal        = std::get<2>(values);
    auto & natural_projection = std::get<3>(values);
    
    this->computeNormalOnElement(element, normal);

    Vector<Real> real_projection(spatial_dimension);
    this->computeProjectionOnElement(element, normal, query,
				     natural_projection, real_projection);

    Vector<Real> distance(spatial_dimension);
    distance = query - real_projection;
    gap = Math::norm(spatial_dimension, distance.storage());

    Vector<Real> direction = distance.normalize(); 
    Real cos_angle = direction.dot(normal);

    Real tolerance = 1e-8;
    
    if (std::abs(cos_angle - 1) <= tolerance && detection_type == _explicit) {
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

  // use contains function to check whether projection lies inside
  // the element, if yes it is a valid projection otherwise no
  // still have to think about what to do if normal exists but
  // projection doesnot lie inside the element 

  bool validity = this->isValidProjection(element, real_projection, natural_projection);

}

/* -------------------------------------------------------------------------- */
bool ContactDetector::isValidProjection(const Element & element, Vector<Real> & real_projection,
					Vector<Real> & natural_projection) {

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type,
					 _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(this->positions /*mesh.getNodes()*/, nodes_coord.storage(),
				     elem_val + element.element * nb_nodes_per_element,
				     nb_nodes_per_element, spatial_dimension);
  
#define GET_NATURAL_COORDINATE(type)					\
  ElementClass<type>::inverseMap(real_projection, nodes_coord, natural_projection) 
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NATURAL_COORDINATE);
#undef GET_NATURAL_COORDINATE

  /*Vector<Real> barycenter(spatial_dimension);
  mesh.getBarycenter(element, barycenter);

  Real distance = 0;
  switch (this->spatial_dimension) {
  case 2: {
    distance = Math::distance_2d(real_projection.storage(), barycenter.storage());
    break;
  }
  case 3: {
    distance = Math::distance_3d(real_projection.storage(), barycenter.storage());    
    break;
  }  
  default: { AKANTU_ERROR("Unknown dimension : " << spatial_dimension); }
  }

  if (distance <= max_dd) {
    return true;
    }*/

  return false;
}

  
/* -------------------------------------------------------------------------- */
void ContactDetector::coordinatesOfElement(const Element & el, Matrix<Real> & coords) {
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
void ContactDetector::normalProjection(const Element & el, const Vector<Real> & slave_coord,
					 Vector<Real> & natural_coord, Real & tolerance) {

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
