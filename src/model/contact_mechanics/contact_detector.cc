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
 ContactDetector::ContactDetector(Mesh & mesh, std::string master, std::string slave, const ID & id,  UInt memory_id)
   : ContactDetector(mesh, mesh.getNodes(), master, slave, id, memory_id) {
  
 }
  
/* -------------------------------------------------------------------------- */
  ContactDetector::ContactDetector(Mesh & mesh, Array<Real> & positions, std::string master, std::string slave, const ID & id, UInt memory_id)
  : Memory(id, memory_id),
    Parsable(ParserType::_contact_detector, id),
    mesh(mesh),
    positions(positions) {

  AKANTU_DEBUG_IN();  

  this->spatial_dimension = mesh.getSpatialDimension();

  this->registerParam("master_surface", master_surface, _pat_parsmod,
		      "Master surface id");
  this->registerParam("slave_surface", slave_surface, _pat_parsmod,
  		      "Slave surface id");

  //auto & ms = this->get("master_surface");

  //const Parser & parser = getStaticParser();
  //const ParserSection & section_detector =
  //    *(parser.getSubSections(ParserType::_contact_detector).first);
  //this->parseSection(section_detector);
  
  //auto sub_sections = getStaticParser();
  //this->master_surface = section_detector.getParameter("master_surface",  _ppsc_current_scope);
  //this->master_surface = section_detector.getParameter("slave_surface",  _ppsc_current_scope);

  this->master_id = master;
  this->slave_id = slave;

  this->mesh.fillNodesToElements();  

  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
void ContactDetector::getMaximalDetectionDistance() {

  AKANTU_DEBUG_IN();

  Real el_size;
  Real max_el_size = std::numeric_limits<Real>::min();

  auto & master_group =
    mesh.getElementGroup(master_surface);

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
void ContactDetector::setMasterSurface(std::string master_surface) {
  this->master_surface = master_surface;  
  this->getMaximalDetectionDistance();
}



/* -------------------------------------------------------------------------- */
void ContactDetector::setSlaveSurface(std::string slave_surface) {
  this->slave_surface = slave_surface;
}

  
/* -------------------------------------------------------------------------- */
void ContactDetector::search() {
  this->globalSearch();
}
  
 
/* -------------------------------------------------------------------------- */
void ContactDetector::globalSearch() {
  
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
  
  AKANTU_DEBUG_INFO( "Intersection BBox "
		     << bbox_intersection );
 
  Vector<Real> center(spatial_dimension);
  bbox_intersection.getCenter(center);

  Vector<Real> spacing(spatial_dimension);
  this->computeCellSpacing(spacing);
    
  auto & master_surface_list =
    mesh.getElementGroup(master_surface).getNodeGroup().getNodes();

  auto & slave_surface_list =
    mesh.getElementGroup(slave_surface).getNodeGroup().getNodes();
 
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
  this->localSearch(slave_grid, master_grid);
}

/* -------------------------------------------------------------------------- */
  void ContactDetector::localSearch(SpatialGrid<UInt> & slave_grid, SpatialGrid<UInt> & master_grid) {

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

   // find the closet master node for each slave node
   for (auto && cell_id : slave_grid) {
     AKANTU_DEBUG_INFO("Looping on next cell");

          
     for (auto && q1: slave_grid.getCell(cell_id)) {
       
       Vector<Real> pos(spatial_dimension);
       for (UInt s: arange(spatial_dimension)) {
	 pos(s) = this->positions(q1, s);
       }

       Real closet_distance = std::numeric_limits<Real>::max();
       UInt closet_master_node;

       // loop over all the neighboring cells of the current cells
       for (auto && neighbor_cell : cell_id.neighbors()) {
	 //AKANTU_DEBUG_INFO("Looping on neighbor cell");
	 // loop over the data of neighboring cells from master grid
	 
	 for (auto && q2 : master_grid.getCell(neighbor_cell)) {
	   
	   AKANTU_DEBUG_INFO("Looping on neighbor cell in master");
	   Vector<Real> pos2(spatial_dimension);
	   for (UInt s: arange(spatial_dimension)) {
	     pos2(s) = this->positions(q2, s);
	   }

	   Real distance = pos.distance(pos2);
	   std::cout << distance << " " << closet_distance << std::endl;
   
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

   
   for (auto & master: master_nodes) {
     Array<Element> elements;
     mesh.getAssociatedElements(master, elements);
     
     //Array<Real> normals(elements.size(), spatial_dimension);
     //this->computeNormalOnElements(elements, normals);

     //Array<Real> projections(elements.size(), 0);
     //this->computeOrthogonalProjections(normals, projections);

     // create a contact element
     // elements array
     // minimum orthogonal projection
     // and the corresponding normal
     // UInt slave node
     // Element master element
     // contact_element = std::make_shared<ContactElement>(UInt, Element);
     // contact_element->setGap(projection)
     // contact_element->setNormal(normal(p));
     // contact_element->setPatch(elements);
     // elements.push_back(contact_element);
     
   }

  // PART III
  // compute the orthogonal distance from each associated element
  // this->computeOrthogonalDistance(*sit, elements);
  // give orthogonal distance
   //   /// compute normals on facets
   // this->computeNormals(); belongs to model class
   // put this function in contact model class
   // and create the mesh_facets in conatct detection class
   // this->getFEEngine("FacetsFEEngine")
   //.computeNormalsOnIntegrationPoints(_not_ghost);
   // const auto & tangents = model->getTangents(type_facet);
  //  const auto & normals = model->getFEEngine("FacetsFEEngine")
  //                             .getNormalsOnIntegrationPoints(type_facet);

   
  // PART IV
  // create contact element,
  // contact element should have one slave node and master facets
  // ContactElement<UInt, Array<Element> > e(slave, elements);
  // set the orthogonal and orthogonal distance
  // 
    
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
      lower_bound(s) -= this->max_dd;
      upper_bound(s) += this->max_dd;
  }
  
  
  AKANTU_DEBUG_INFO("BBox" << bbox);
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeCellSpacing(Vector<Real> & spacing) {

  for (UInt s: arange(spatial_dimension)) 
    spacing(s) = std::sqrt(2.0) * max_dd;
    
}
  
/* -------------------------------------------------------------------------- */
void ContactDetector::computeNormalOnElements(Array<Element> & elements, Array<Real> & normals) {

  // get from connectivity the nodes for each element in Array elements;

  const Matrix<Real> & coords;
  this->coordinatesOfElement(e, coords);
  
  Matrix J(this->spatial_dimension, this->spatial_dimension);
  this->vectorsAlongElement(coords, J);  

  switch (this->spatial_dimension) {
  case 2: {
    // extract the vectors
    Math::normal2(J.storage(), normal.storage());
    break;
  }
  case 3: {
    Math::normal3(J(0).storage(), J(1).storage, normal(p).storage());
    break;
  }  
  default:
    break;
  }
}

  

/* -------------------------------------------------------------------------- */
/* To be put in conatct mechancis model class

void ContactDetector::computeNormals() {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = this->inserter->getMeshFacets();
  this->getFEEngine("FacetsFEEngine")
      .computeNormalsOnIntegrationPoints(_not_ghost);

  UInt tangent_components =
      Model::spatial_dimension * (Model::spatial_dimension - 1);

  tangents.initialize(mesh_facets, _nb_component = tangent_components,
                      _spatial_dimension = Model::spatial_dimension - 1);
 
  for (auto facet_type :
       mesh_facets.elementTypes(Model::spatial_dimension - 1)) {
    const Array<Real> & normals =
        this->getFEEngine("FacetsFEEngine")
            .getNormalsOnIntegrationPoints(facet_type);

    Array<Real> & tangents = this->tangents(facet_type);

    Math::compute_tangents(normals, tangents);
  }

  AKANTU_DEBUG_OUT();
}*/

  
  
} // akantu
