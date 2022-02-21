/**
 * @file   contact_detector.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Dec 05 2018
 * @date last modification: Thu Jun 24 2021
 *
 * @brief  Mother class for all detection algorithms
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_detector.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactDetector::ContactDetector(Mesh & mesh, const ID & id)
    : ContactDetector(mesh, mesh.getNodes(), id) {}

/* -------------------------------------------------------------------------- */
ContactDetector::ContactDetector(Mesh & mesh, Array<Real> positions,
                                 const ID & id)
    : Parsable(ParserType::_contact_detector, id), mesh(mesh),
      positions(0, mesh.getSpatialDimension()) {

  AKANTU_DEBUG_IN();

  this->spatial_dimension = mesh.getSpatialDimension();

  this->positions.copy(positions);

  const Parser & parser = getStaticParser();
  const ParserSection & section =
      *(parser.getSubSections(ParserType::_contact_detector).first);

  this->parseSection(section);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactDetector::parseSection(const ParserSection & section) {
  auto type = section.getParameterValue<std::string>("type");

  if (type == "implicit") {
    this->detection_type = _implicit;
  } else if (type == "explicit") {
    this->detection_type = _explicit;
  } else {
    AKANTU_ERROR("Unknown detection type : " << type);
  }

  this->projection_tolerance =
      section.getParameterValue<Real>("projection_tolerance");
  this->max_iterations = section.getParameterValue<Real>("max_iterations");
  this->extension_tolerance =
      section.getParameterValue<Real>("extension_tolerance");
}

/* -------------------------------------------------------------------------- */
void ContactDetector::search(Array<ContactElement> & elements,
                             Array<Real> & gaps, Array<Real> & normals,
                             Array<Real> & tangents,
                             Array<Real> & projections) {
  auto surface_dimension = spatial_dimension - 1;

  this->mesh.fillNodesToElements(surface_dimension);
  this->computeMaximalDetectionDistance();

  contact_pairs.clear();

  SpatialGrid<UInt> master_grid(spatial_dimension);
  SpatialGrid<UInt> slave_grid(spatial_dimension);

  this->globalSearch(slave_grid, master_grid);

  this->localSearch(slave_grid, master_grid);

  this->createContactElements(elements, gaps, normals, tangents, projections);
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
      for (UInt s : arange(spatial_dimension)) {
        pos(s) = this->positions(slave_node, s);
      }

      Real closet_distance = std::numeric_limits<Real>::max();
      UInt closet_master_node;

      /// loop over all the neighboring cells of the current cell
      for (auto && neighbor_cell : cell_id.neighbors()) {
        /// loop over the data of neighboring cells from master grid
        for (auto && master_node : master_grid.getCell(neighbor_cell)) {

          /// check for self contact
          if (slave_node == master_node) {
            continue;
          }

          bool is_valid = true;
          Array<Element> elements;
          this->mesh.getAssociatedElements(slave_node, elements);

          for (auto & elem : elements) {
            if (elem.kind() != _ek_regular) {
              continue;
            }

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
          for (UInt s : arange(spatial_dimension)) {
            pos2(s) = this->positions(master_node, s);
          }

          Real distance = pos.distance(pos2);

          if (distance <= closet_distance and is_valid) {
            closet_master_node = master_node;
            closet_distance = distance;
            pair_exists = true;
          }
        }
      }

      if (pair_exists) {
        contact_pairs.emplace_back(
            std::make_pair(slave_node, closet_master_node));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void ContactDetector::createContactElements(
    Array<ContactElement> & contact_elements, Array<Real> & gaps,
    Array<Real> & normals, Array<Real> & tangents, Array<Real> & projections) {
  auto surface_dimension = spatial_dimension - 1;

  Real alpha;
  switch (detection_type) {
  case _explicit: {
    alpha = 1.0;
    break;
  }
  case _implicit: {
    alpha = -1.0;
    break;
  }
  default:
    AKANTU_EXCEPTION(detection_type
                     << " is not a valid contact detection type");
    break;
  }

  for (auto & pairs : contact_pairs) {

    const auto & slave_node = pairs.first;

    Vector<Real> slave(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      slave(s) = this->positions(slave_node, s);
    }

    const auto & master_node = pairs.second;
    Array<Element> elements;
    this->mesh.getAssociatedElements(master_node, elements);

    auto & gap = gaps.begin()[slave_node];
    Vector<Real> normal(normals.begin(spatial_dimension)[slave_node]);
    Vector<Real> projection(projections.begin(surface_dimension)[slave_node]);
    Matrix<Real> tangent(
        tangents.begin(surface_dimension, spatial_dimension)[slave_node]);
    auto index = GeometryUtils::orthogonalProjection(
        mesh, positions, slave, elements, gap, projection, normal, tangent,
        alpha, this->max_iterations, this->projection_tolerance,
        this->extension_tolerance);

    // if a valid projection is not found on the patch of elements
    // index is -1 or if not a valid self contact, the contact element
    // is not created
    if (index == UInt(-1) or !isValidSelfContact(slave_node, gap, normal)) {
      gap *= 0;
      normal *= 0;
      projection *= 0;
      tangent *= 0;
      continue;
    }

    // create contact element
    contact_elements.push_back(ContactElement(slave_node, elements[index]));
  }

  contact_pairs.clear();
}

} // namespace akantu
