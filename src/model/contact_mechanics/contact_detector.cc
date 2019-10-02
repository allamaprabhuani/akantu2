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
ContactDetector::ContactDetector(Mesh & mesh, const ID & id,
                                 const MemoryID & memory_id)
    : ContactDetector(mesh, mesh.getNodes(), id, memory_id) {}

/* -------------------------------------------------------------------------- */
ContactDetector::ContactDetector(Mesh & mesh, Array<Real> positions,
                                 const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), Parsable(ParserType::_contact_detector, id),
      mesh(mesh), positions(0, mesh.getSpatialDimension()) {

  AKANTU_DEBUG_IN();

  this->spatial_dimension = mesh.getSpatialDimension();

  this->positions.copy(positions);

  this->parseSection();

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
  } else if (type == "explicit") {
    this->detection_type = _explicit;
  } else {
    AKANTU_ERROR("Unknown detection type : " << type);
  }
}

/* -------------------------------------------------------------------------- */
void ContactDetector::search(std::map<UInt, ContactElement> & contact_map) {

  this->mesh.fillNodesToElements(this->spatial_dimension - 1);
  this->computeMaximalDetectionDistance();

  contact_pairs.clear();

  SpatialGrid<UInt> master_grid(spatial_dimension);
  SpatialGrid<UInt> slave_grid(spatial_dimension);

  this->globalSearch(slave_grid, master_grid);

  this->localSearch(slave_grid, master_grid);

  this->constructContactMap(contact_map);
}

/* -------------------------------------------------------------------------- */
void ContactDetector::globalSearch(SpatialGrid<UInt> & slave_grid,
                                   SpatialGrid<UInt> & master_grid) {

  auto & master_list = surface_selector->getMasterList();
  auto & slave_list = surface_selector->getSlaveList();

  BBox bbox_master(spatial_dimension);
  this->constructBoundingBox(bbox_master, master_list);

  BBox bbox_slave(spatial_dimension);
  this->constructBoundingBox(bbox_slave, slave_list);

  auto && bbox_intersection = bbox_master.intersection(bbox_slave);

  AKANTU_DEBUG_INFO("Intersection BBox " << bbox_intersection);

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

  /// find the closet master node for each slave node
  for (auto && cell_id : slave_grid) {
    /// loop over all the slave nodes of the current cell
    for (auto && slave_node : slave_grid.getCell(cell_id)) {

      bool pair_exists = false;

      Vector<Real> pos(spatial_dimension);
      for (UInt s : arange(spatial_dimension))
        pos(s) = this->positions(slave_node, s);

      Real closet_distance = std::numeric_limits<Real>::max();
      UInt closet_master_node;

      /// loop over all the neighboring cells of the current cell
      for (auto && neighbor_cell : cell_id.neighbors()) {
        /// loop over the data of neighboring cells from master grid
        for (auto && master_node : master_grid.getCell(neighbor_cell)) {

          /// check for self contact
          if (slave_node == master_node)
            continue;

	  bool is_valid = true;
          Array<Element> elements;
          this->mesh.getAssociatedElements(slave_node, elements);

          for (auto & elem : elements) {
            if (elem.kind() != _ek_regular)
              continue;

            Vector<UInt> connectivity =
                const_cast<const Mesh &>(this->mesh).getConnectivity(elem);

            auto node_iter = std::find(connectivity.begin(), connectivity.end(),
                                       master_node);
            if (node_iter != connectivity.end()) {
              is_valid = false;
	      break;
            }
          }

          Vector<Real> pos2(spatial_dimension);
          for (UInt s : arange(spatial_dimension))
            pos2(s) = this->positions(master_node, s);

          Real distance = pos.distance(pos2);

          if (distance <= closet_distance and is_valid) {
            closet_master_node = master_node;
            closet_distance = distance;
            pair_exists = true;
          }
        }
      }

      if (pair_exists)
        contact_pairs.push_back(std::make_pair(slave_node, closet_master_node));
    }
  }
}

/* -------------------------------------------------------------------------- */
void ContactDetector::constructContactMap(
    std::map<UInt, ContactElement> & contact_map) {

  auto surface_dimension = spatial_dimension - 1;

  std::map<UInt, ContactElement> previous_contact_map = contact_map;
  contact_map.clear();

  auto get_connectivity = [&](auto & slave, auto & master) {
    Vector<UInt> master_conn =
        const_cast<const Mesh &>(this->mesh).getConnectivity(master);
    Vector<UInt> elem_conn(master_conn.size() + 1);

    elem_conn[0] = slave;
    for (UInt i = 1; i < elem_conn.size(); ++i) {
      elem_conn[i] = master_conn[i - 1];
    }

    return elem_conn;
  };

  for (auto & pairs : contact_pairs) {

    const auto & slave_node = pairs.first;
    const auto & master_node = pairs.second;

    Array<Element> all_elements;
    this->mesh.getAssociatedElements(master_node, all_elements);

    Array<Element> boundary_elements;
    this->filterBoundaryElements(all_elements, boundary_elements);

    Array<Real> gaps(boundary_elements.size(), 1, "gaps");
    Array<Real> normals(boundary_elements.size(), spatial_dimension, "normals");
    Array<Real> projections(boundary_elements.size(), surface_dimension,
                            "projections");

    auto index = this->computeOrthogonalProjection(
        slave_node, boundary_elements, normals, gaps, projections);

    if (index == UInt(-1)) {
      continue;
    }

    auto connectivity = get_connectivity(slave_node, boundary_elements[index]);

    // assign contact element attributes
    contact_map[slave_node].setMaster(boundary_elements[index]);
    contact_map[slave_node].setGap(gaps[index]);
    contact_map[slave_node].setNormal(
        Vector<Real>(normals.begin(spatial_dimension)[index], true));
    contact_map[slave_node].setProjection(
        Vector<Real>(projections.begin(surface_dimension)[index], true));
    contact_map[slave_node].setConnectivity(connectivity);

    // tangent computation on master surface
    Matrix<Real> tangents(surface_dimension, spatial_dimension);
    this->computeTangentsOnElement(contact_map[slave_node].master,
                                   contact_map[slave_node].projection,
                                   contact_map[slave_node].normal, tangents);

    contact_map[slave_node].setTangent(tangents);

    // friction calculation requires history of previous natural
    // projection as well as traction
    auto search = previous_contact_map.find(slave_node);
    if (search != previous_contact_map.end()) {
      auto previous_projection =
          previous_contact_map[slave_node].getPreviousProjection();
      contact_map[slave_node].setPreviousProjection(previous_projection);

      auto previous_traction = previous_contact_map[slave_node].getTraction();
      contact_map[slave_node].setTraction(previous_traction);
    } else {
      Vector<Real> previous_projection(surface_dimension, 0.);
      contact_map[slave_node].setPreviousProjection(previous_projection);

      Vector<Real> previous_traction(surface_dimension, 0.);
      contact_map[slave_node].setTraction(previous_traction);
    }

    // to ensure the self contact between surface does not lead to
    // detection of master element which is closet but not
    // orthogonally opposite to the slave surface
    bool is_valid_self_contact =
        this->checkValidityOfSelfContact(slave_node, contact_map[slave_node]);

    if (!is_valid_self_contact) {
      contact_map.erase(slave_node);
    }
  }

  contact_pairs.clear();
}

/* -------------------------------------------------------------------------- */
UInt ContactDetector::computeOrthogonalProjection(
    const UInt & node, const Array<Element> & elements, Array<Real> & normals,
    Array<Real> & gaps, Array<Real> & projections) {

  Vector<Real> query(spatial_dimension);
  for (UInt s : arange(spatial_dimension))
    query(s) = this->positions(node, s);

  UInt counter = 0;
  UInt index = UInt(-1);

  Real min_gap = std::numeric_limits<Real>::max();

  for (auto && values :
       zip(elements, gaps, make_view(normals, spatial_dimension),
           make_view(projections, spatial_dimension - 1))) {

    const auto & element = std::get<0>(values);
    auto & gap = std::get<1>(values);
    auto & normal = std::get<2>(values);
    auto & projection = std::get<3>(values);

    this->computeNormalOnElement(element, normal);

    Vector<Real> real_projection(spatial_dimension);
    this->computeProjectionOnElement(element, normal, query, projection,
                                     real_projection);

    gap = this->computeGap(query, real_projection);

    // check if gap is valid or not
    // to check this we need normal on master element, vector from
    // real projection to slave node, if it is explicit detection
    // scheme, we want it to interpenetrate, the dot product should be
    // negative and -1.0 and for implciit detection , the dot product
    // should be positive and 1.0

    bool is_valid = this->checkValidityOfProjection(projection);

    auto master_to_slave = query - real_projection;
    auto norm = master_to_slave.norm();
    if (norm != 0)
      master_to_slave /= norm;

    Real tolerance = 1e-8;

    switch (detection_type) {
    case _explicit: {
      auto product = master_to_slave.dot(normal);
      auto variation = std::abs(product + 1.0);

      if (variation <= tolerance and gap <= min_gap and is_valid) {
        min_gap = gap;
        index = counter;
        gap *= -1.0;
      }
      break;
    }
    case _implicit: {
      auto product = master_to_slave.dot(normal);
      auto variation = std::abs(product - 1.0);

      if (variation <= tolerance and gap <= min_gap and is_valid) {
        min_gap = gap;
        index = counter;
        gap *= -1.0;
      }
    }
    default:
      break;
    }

    counter++;
  }

  return index;
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeProjectionOnElement(
    const Element & element, const Vector<Real> & normal,
    const Vector<Real> & query, Vector<Real> & natural_projection,
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
void ContactDetector::computeNaturalProjection(
    const Element & element, Vector<Real> & real_projection,
    Vector<Real> & natural_projection) {

  const ElementType & type = element.type;
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt * elem_val = mesh.getConnectivity(type, _not_ghost).storage();
  Matrix<Real> nodes_coord(spatial_dimension, nb_nodes_per_element);

  mesh.extractNodalValuesFromElement(this->positions, nodes_coord.storage(),
                                     elem_val +
                                         element.element * nb_nodes_per_element,
                                     nb_nodes_per_element, spatial_dimension);

#define GET_NATURAL_COORDINATE(type)                                           \
  ElementClass<type>::inverseMap(real_projection, nodes_coord,                 \
                                 natural_projection)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_NATURAL_COORDINATE);
#undef GET_NATURAL_COORDINATE
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeTangentsOnElement(const Element & el,
                                               Vector<Real> & projection,
                                               Matrix<Real> & tangents) {

  const ElementType & type = el.type;

  UInt nb_nodes_master = Mesh::getNbNodesPerElement(type);

  Vector<Real> shapes(nb_nodes_master);
  Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);

#define GET_SHAPES_NATURAL(type)                                               \
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

#define GET_SHAPE_DERIVATIVES_NATURAL(type)                                    \
  ElementClass<type>::computeDNDS(projection, shapes_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

  Matrix<Real> coords(spatial_dimension, nb_nodes_master);
  coordinatesOfElement(el, coords);

  tangents.mul<false, true>(shapes_derivatives, coords);

  auto temp_tangents = tangents.transpose();
  for (UInt i = 0; i < spatial_dimension - 1; ++i) {
    auto temp = Vector<Real>(temp_tangents(i));
    temp_tangents(i) = temp.normalize();
  }

  tangents = temp_tangents.transpose();
}

/* -------------------------------------------------------------------------- */
void ContactDetector::computeTangentsOnElement(const Element & el,
                                               Vector<Real> & projection,
                                               Vector<Real> & normal,
                                               Matrix<Real> & tangents) {

  const ElementType & type = el.type;

  UInt nb_nodes_master = Mesh::getNbNodesPerElement(type);

  Vector<Real> shapes(nb_nodes_master);
  Matrix<Real> shapes_derivatives(spatial_dimension - 1, nb_nodes_master);

#define GET_SHAPES_NATURAL(type)                                               \
  ElementClass<type>::computeShapes(projection, shapes)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPES_NATURAL);
#undef GET_SHAPES_NATURAL

#define GET_SHAPE_DERIVATIVES_NATURAL(type)                                    \
  ElementClass<type>::computeDNDS(projection, shapes_derivatives)
  AKANTU_BOOST_ALL_ELEMENT_SWITCH(GET_SHAPE_DERIVATIVES_NATURAL);
#undef GET_SHAPE_DERIVATIVES_NATURAL

  Matrix<Real> coords(spatial_dimension, nb_nodes_master);
  coordinatesOfElement(el, coords);

  tangents.mul<false, true>(shapes_derivatives, coords);

  auto temp_tangents = tangents.transpose();
  for (UInt i = 0; i < spatial_dimension - 1; ++i) {
    auto temp = Vector<Real>(temp_tangents(i));
    temp_tangents(i) = temp.normalize();
  }

  tangents = temp_tangents.transpose();

  // to ensure that direction of tangents are correct, cross product
  // of tangents should give the normal vector computed earlier
  switch (spatial_dimension) {
  case 2: {
    Vector<Real> e_z(3);
    e_z[0] = 0.;
    e_z[1] = 0.;
    e_z[2] = 1.;

    Vector<Real> tangent(3);
    tangent[0] = tangents(0, 0);
    tangent[1] = tangents(0, 1);
    tangent[2] = 0.;

    auto exp_normal = e_z.crossProduct(tangent);

    auto & cal_normal = normal;

    auto ddot = cal_normal.dot(exp_normal);
    if (ddot < 0) {
      tangents *= -1.0;
    }

    break;
  }
  case 3: {
    auto tang_trans = tangents.transpose();
    auto tang1 = Vector<Real>(tang_trans(0));
    auto tang2 = Vector<Real>(tang_trans(1));

    auto tang1_cross_tang2 = tang1.crossProduct(tang2);
    auto exp_normal = tang1_cross_tang2 / tang1_cross_tang2.norm();

    auto & cal_normal = normal;

    auto ddot = cal_normal.dot(exp_normal);
    if (ddot < 0) {
      tang_trans(1) *= -1.0;
    }

    tangents = tang_trans.transpose();
    break;
  }
  default:
    break;
  }
}

} // namespace akantu
