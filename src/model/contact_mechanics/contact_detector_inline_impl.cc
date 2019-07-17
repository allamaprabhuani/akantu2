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
inline UInt ContactDetector::getElementIndex(Array<Real> & gaps,
                                             Array<Real> & projections,
					     Array<Real> & normals) {
  auto surface_dimension = spatial_dimension - 1;

  UInt index;
  Real gap_min = std::numeric_limits<Real>::max();

  UInt counter = 0;
  for (auto && values : zip(gaps, make_view(projections, surface_dimension))) {
    auto & gap = std::get<0>(values);
    auto & projection = std::get<1>(values);

    bool is_valid = this->checkValidityOfProjection(projection);

    if (is_valid and gap <= gap_min) {
      gap_min = gap;
      index = counter;
    }
    counter++;
  }

  // if a slave node does not have a valid projection on any of the
  // element, the minimum projection is taken as valid projection,
  // this sometimes leads to situation where a slave node is outside
  // but is still consider for potential contact
  if (index >= gaps.size()) {
    auto gap_min_it = std::min_element(gaps.begin(), gaps.end());
    auto index_it = std::find(gaps.begin(), gaps.end(), *gap_min_it);
    index = *index_it;
  }

  return index;
}

/* -------------------------------------------------------------------------- */
inline bool
ContactDetector::checkValidityOfProjection(Vector<Real> & projection) {

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
inline void ContactDetector::coordinatesOfElement(const Element & el,
                                                  Matrix<Real> & coords) {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(el.type);
  Vector<UInt> connect = mesh.getConnectivity(el.type, _not_ghost)
                             .begin(nb_nodes_per_element)[el.element];

  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt node = connect[n];
    for (UInt s : arange(spatial_dimension)) {
      coords(s, n) = this->positions(node, s);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeCellSpacing(Vector<Real> & spacing) {
  for (UInt s : arange(spatial_dimension))
    spacing(s) = std::sqrt(2.0) * max_dd;
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::constructBoundingBox(BBox & bbox,
                                      const Array<UInt> & nodes_list) {

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

  for (UInt s : arange(spatial_dimension)) {
    lower_bound(s) -= this->max_bb;
    upper_bound(s) += this->max_bb;
  }

  AKANTU_DEBUG_INFO("BBox" << bbox);
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::constructGrid(SpatialGrid<UInt> & grid,
                                           BBox & bbox,
                                           const Array<UInt> & nodes_list) {
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
template <Surface id> inline std::string ContactDetector::getSurfaceId() {
  return surfaces[id];
}

/* -------------------------------------------------------------------------- */
template <Surface id>
inline void ContactDetector::setSurfaceId(const std::string name) {
  surfaces[id] = name;
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeMaximalDetectionDistance() {

  AKANTU_DEBUG_IN();

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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline Vector<UInt>
ContactDetector::constructConnectivity(UInt & slave, const Element & master) {

  Vector<UInt> master_conn =
      const_cast<const Mesh &>(this->mesh).getConnectivity(master);

  Vector<UInt> elem_conn(master_conn.size() + 1);

  elem_conn[0] = slave;
  for (UInt i = 1; i < elem_conn.size(); ++i) {
    elem_conn[i] = master_conn[i - 1];
  }

  return elem_conn;
}

/* -------------------------------------------------------------------------- */
inline void ContactDetector::computeNormalOnElement(const Element & element,
                                                    Vector<Real> & normal) {

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
                                                 Matrix<Real> & vectors) {

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
inline Real ContactDetector::computeGap(Vector<Real> & slave,
                                        Vector<Real> & master,
                                        Vector<Real> & normal) {

  Vector<Real> slave_to_master(spatial_dimension);
  slave_to_master = master - slave;

  Real gap = slave_to_master.norm();

  auto projection = slave_to_master.dot(normal);

  // if slave node is beneath the master surface
  if (projection > 0) {
    //gap *= -1.0;
  }

  return gap;
}

/* -------------------------------------------------------------------------- */
inline void
ContactDetector::filterBoundaryElements(Array<Element> & elements,
                                        Array<Element> & boundary_elements) {

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
      if (subelem == ElementNull)
        continue;

      if (subelem.kind() == _ek_regular)
        ++nb_subelements_regular;
    }

    auto nb_subelements = element_to_subelement.size();

    if (nb_subelements_regular < nb_subelements)
      boundary_elements.push_back(elem);
  }
}

/* -------------------------------------------------------------------------- */
inline bool
ContactDetector::checkValidityOfSelfContact(const UInt & slave_node,
                                            const ContactElement & element) {

  UInt master_node;

  for (auto & pair : contact_pairs) {
    if (pair.first == slave_node) {
      master_node = pair.second;
      break;
    }
  }

  Array<Element> elements;
  this->mesh.getAssociatedElements(slave_node, elements);
  
  for (auto & elem : elements) {
    if (elem.kind() != _ek_regular)
      continue;

    Vector<UInt> connectivity =
        const_cast<const Mesh &>(this->mesh).getConnectivity(elem);

    auto node_iter =
        std::find(connectivity.begin(), connectivity.end(), master_node);
    if (node_iter != connectivity.end()) {
      return false;
    }
  }

  if (std::abs(element.gap) > 2.0 * min_dd) {
    return false;
  }

  return true;
}

} // namespace akantu

#endif /*  __AKANTU_CONTACT_DETECTOR_INLINE_IMPL_CC__ */
