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
template <class Derived, std::enable_if_t<aka::is_vector<Derived>::value> *>
inline bool ContactDetector::checkValidityOfProjection(
    Eigen::MatrixBase<Derived> & projection) const {
  Real tolerance = 1e-3;
  return std::all_of(projection.begin(), projection.end(),
                     [&tolerance](auto && xi) {
                       return (xi > -1.0 - tolerance) or (xi < 1.0 + tolerance);
                     });
}

/* -------------------------------------------------------------------------- */
template <class Derived>
inline void ContactDetector::coordinatesOfElement(
    const Element & el, Eigen::MatrixBase<Derived> & coords) const {
  coords = mesh.extractNodalValuesFromElement(positions, el);
}

/* -------------------------------------------------------------------------- */
template <class Derived, std::enable_if_t<aka::is_vector<Derived>::value> *>
inline void ContactDetector::computeCellSpacing(
    Eigen::MatrixBase<Derived> & spacing) const {
  spacing.fill(std::sqrt(2.0) * max_dd);
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::constructGrid(SpatialGrid<Idx> & grid, BBox & bbox,
                               const Array<Idx> & nodes_list) const {
  auto position_it = make_view(this->positions, spatial_dimension).begin();
  auto to_grid = [&](auto node) {
    auto && pos = position_it[node];
    if (bbox.contains(pos)) {
      grid.insert(node, pos);
    }
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_grid);
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::constructBoundingBox(BBox & bbox,
                                      const Array<Idx> & nodes_list) const {
  auto to_bbox = [&](auto node) {
    auto && pos = make_view(this->positions, spatial_dimension).begin()[node];
    bbox += pos;
  };

  std::for_each(nodes_list.begin(), nodes_list.end(), to_bbox);

  Vector<Real> lower_bound = bbox.getLowerBounds().array() - max_bb;
  Vector<Real> upper_bound = bbox.getUpperBounds().array() + max_bb;

  bbox.setLowerBounds(lower_bound);
  bbox.setUpperBounds(upper_bound);

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
      Int nb_nodes_per_element = mesh.getNbNodesPerElement(element.type);
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
inline Vector<Idx>
ContactDetector::constructConnectivity(Idx & slave,
                                       const Element & master) const {
  auto && master_conn = this->mesh.getConnectivity(master);

  Vector<Idx> elem_conn(master_conn.size() + 1);
  elem_conn[0] = slave;
  elem_conn.block(1, 0, master_conn.size(), 1) = master_conn;

  return elem_conn;
}

/* -------------------------------------------------------------------------- */
template <class Derived, std::enable_if_t<aka::is_vector<Derived>::value> *>
inline void ContactDetector::computeNormalOnElement(
    const Element & element, Eigen::MatrixBase<Derived> & normal) const {
  Matrix<Real> vectors(spatial_dimension, spatial_dimension - 1);
  this->vectorsAlongElement(element, vectors);

  switch (this->spatial_dimension) {
  case 2: {
    normal = Math::normal(vectors);
    break;
  }
  case 3: {
    normal = Math::normal(vectors.col(0), vectors.col(1));
    break;
  }
  default: {
    AKANTU_ERROR("Unknown dimension : " << spatial_dimension);
  }
  }

  // to ensure that normal is always outwards from master surface
  const auto & element_to_subelement =
      mesh.getElementToSubelement(element.type)(element.element);

  Vector<Real> outside = mesh.getBarycenter(element);

  // check if mesh facets exists for cohesive elements contact
  Vector<Real> inside;
  if (mesh.isMeshFacets()) {
    inside = mesh.getMeshParent().getBarycenter(element_to_subelement[0]);
  } else {
    inside = mesh.getBarycenter(element_to_subelement[0]);
  }

  auto projection = (outside - inside).dot(normal);

  if (projection < 0) {
    normal *= -1.0;
  }
}

/* -------------------------------------------------------------------------- */
template <class Derived>
inline void ContactDetector::vectorsAlongElement(
    const Element & el, Eigen::MatrixBase<Derived> & vectors) const {
  auto nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);

  Matrix<Real> coords(spatial_dimension, nb_nodes_per_element);
  this->coordinatesOfElement(el, coords);

  for (auto i : arange(spatial_dimension - 1)) {
    vectors(i) = coords(i + 1) - coords(0);
  }
}

/* -------------------------------------------------------------------------- */
template <class Derived1, class Derived2,
          std::enable_if_t<aka::are_vectors<Derived1, Derived2>::value> *>
inline Real
ContactDetector::computeGap(const Eigen::MatrixBase<Derived1> & slave,
                            const Eigen::MatrixBase<Derived2> & master) const {
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
template <class Derived, std::enable_if_t<aka::is_vector<Derived>::value> *>
inline bool ContactDetector::isValidSelfContact(
    const Idx & slave_node, const Real & gap,
    const Eigen::MatrixBase<Derived> & normal) const {
  Idx master_node{-1};

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

    auto && connectivity = this->mesh.getConnectivity(element);

    auto coords = mesh.extractNodalValuesFromElement(positions, element);

    // finding the normal at slave node by averaging of normals
    auto normal = GeometryUtils::normal(mesh, coords, element);
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
