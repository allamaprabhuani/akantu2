/**
 * @file   contact_detector_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed May 08 2019
 * @date last modification: Thu Jun 24 2021
 *
 * @brief  inine implementation of the contact detector class
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

#ifndef __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__
#define __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline bool
ContactDetector::checkValidityOfProjection(Vector<Real> & projection) const {
  Real tolerance = 1e-3;
  return std::all_of(projection.begin(), projection.end(),
                     [&tolerance](auto && xi) {
                       return (xi > -1.0 - tolerance) or (xi < 1.0 + tolerance);
                     });
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::coordinatesOfElement(const Element & el,
                                                  Matrix<Real> & coords) const {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  const Vector<UInt> connect = mesh.getConnectivity(el.type, _not_ghost)
                                   .begin(nb_nodes_per_element)[el.element];

  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s : arange(spatial_dimension)) {
      coords(s, n) = this->positions(node, s);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeCellSpacing(Vector<Real> & spacing) const {
  for (UInt s : arange(spatial_dimension)) {
    spacing(s) = std::sqrt(2.0) * max_dd;
  }
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::constructGrid(SpatialGrid<UInt> & grid, BBox & bbox,
                               const Array<UInt> & nodes_list) const {
  auto to_grid = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      pos(s) = this->positions(node, s);
    }

    if (bbox.contains(pos)) {
      grid.insert(node, pos);
    }
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_grid);
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::constructBoundingBox(BBox & bbox,
                                      const Array<UInt> & nodes_list) const {
  auto to_bbox = [&](UInt node) {
    Vector<Real> pos(spatial_dimension);
    for (UInt s : arange(spatial_dimension)) {
      pos(s) = this->positions(node, s);
    }
    bbox += pos;
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_bbox);

  auto & lower_bound = bbox.getLowerBounds();
  auto & upper_bound = bbox.getUpperBounds();

  lower_bound -= this->max_bb;
  upper_bound += this->max_bb;

  AKANTU_DEBUG_INFO("BBox" << bbox);
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeMaximalDetectionDistance() {
  Real elem_size;
  Real max_elem_size = std::numeric_limits<Real>::min();
  Real min_elem_size = std::numeric_limits<Real>::max();

  auto & master_nodes = this->surface_selector->getMasterList();

  for (auto & master : master_nodes) {
    Array<Element> elements;
    this->mesh.getAssociatedElements(master, elements);

    for (auto element : elements) {
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element.type);
      Matrix<Real> elem_coords(spatial_dimension, nb_nodes_per_element);
      this->coordinatesOfElement(element, elem_coords);

      elem_size = FEEngine::getElementInradius(elem_coords, element.type);
      max_elem_size = std::max(max_elem_size, elem_size);
      min_elem_size = std::min(min_elem_size, elem_size);
    }
  }

  AKANTU_DEBUG_INFO("The maximum element size : " << max_elem_size);

  this->min_dd = min_elem_size;
  this->max_dd = max_elem_size;
  this->max_bb = max_elem_size;
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt>
ContactDetector::constructConnectivity(UInt & slave,
                                       const Element & master) const {
  const Vector<UInt> master_conn = this->mesh.getConnectivity(master);

  Vector<UInt> elem_conn(master_conn.size() + 1);
  elem_conn[0] = slave;
  for (UInt i = 1; i < elem_conn.size(); ++i) {
    elem_conn[i] = master_conn[i - 1];
  }

  return elem_conn;
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::computeNormalOnElement(const Element & element,
                                        Vector<Real> & normal) const {
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
  default: {
    AKANTU_ERROR("Unknown dimension : " << spatial_dimension);
  }
  }

  // to ensure that normal is always outwards from master surface
  const auto & element_to_subelement =
      mesh.getElementToSubelement(element.type)(element.element);

  Vector<Real> outside(spatial_dimension);
  mesh.getBarycenter(element, outside);

  // check if mesh facets exists for cohesive elements contact
  Vector<Real> inside(spatial_dimension);
  if (mesh.isMeshFacets()) {
    mesh.getMeshParent().getBarycenter(element_to_subelement[0], inside);
  } else {
    mesh.getBarycenter(element_to_subelement[0], inside);
  }

  Vector<Real> inside_to_outside = outside - inside;
  auto projection = inside_to_outside.dot(normal);

  if (projection < 0) {
    normal *= -1.0;
  }
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::vectorsAlongElement(const Element & el,
                                                 Matrix<Real> & vectors) const {
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  this->coordinatesOfElement(el, coords);

  for (auto i : arange(spatial_dimension - 1)) {
    vectors(i) = Vector<Real>(coords(i + 1)) - Vector<Real>(coords(0));
  }
}

/* -------------------------------------------------------------------------- */
inline Real ContactDetector::computeGap(const Vector<Real> & slave,
                                        const Vector<Real> & master) const {
  auto gap = (master - slave).norm();
  return gap;
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::filterBoundaryElements(
    const Array<Element> & elements, Array<Element> & boundary_elements) const {
  for (auto elem : elements) {
    const auto & element_to_subelement =
        mesh.getElementToSubelement(elem.type)(elem.element);

    // for regular boundary elements
    if (element_to_subelement.size() == 1 and
        element_to_subelement[0].kind() == _ek_regular) {
      boundary_elements.push_back(elem);
      continue;
    }

    // for cohesive boundary elements
    UInt nb_subelements_regular = 0;
    for (auto subelem : element_to_subelement) {
      if (subelem == ElementNull) {
        continue;
      }

      if (subelem.kind() == _ek_regular) {
        ++nb_subelements_regular;
      }
    }

    auto nb_subelements = element_to_subelement.size();

    if (nb_subelements_regular < nb_subelements) {
      boundary_elements.push_back(elem);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline bool
ContactDetector::isValidSelfContact(const UInt & slave_node, const Real & gap,
                                    const Vector<Real> & normal) const {
  UInt master_node;

  // finding the master node corresponding to slave node
  for (auto && pair : contact_pairs) {
    if (pair.first == slave_node) {
      master_node = pair.second;
      break;
    }
  }

  Array<Element> slave_elements;
  this->mesh.getAssociatedElements(slave_node, slave_elements);

  // Check 1 : master node is not connected to elements connected to
  // slave node
  Vector<Real> slave_normal(spatial_dimension);
  for (auto & element : slave_elements) {
    if (element.kind() != _ek_regular) {
      continue;
    }

    const Vector<UInt> connectivity = this->mesh.getConnectivity(element);

    // finding the normal at slave node by averaging of normals
    Vector<Real> normal(spatial_dimension);
    GeometryUtils::normal(mesh, positions, element, normal);
    slave_normal = slave_normal + normal;

    auto node_iter =
        std::find(connectivity.begin(), connectivity.end(), master_node);
    if (node_iter != connectivity.end()) {
      return false;
    }
  }

  // Check 2 : if gap is twice the size of smallest element
  if (std::abs(gap) > 2.0 * min_dd) {
    return false;
  }

  // Check 3 : check the directions of normal at slave node and at
  // master element, should be in opposite directions
  auto norm = slave_normal.norm();
  if (norm != 0) {
    slave_normal /= norm;
  }

  auto product = slave_normal.dot(normal);

  return not(product >= 0);
}

} // namespace akantu

#endif /*  __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__ */
