/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
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
 */

/* -------------------------------------------------------------------------- */
#include "neighborhood_max_criterion.hh"
#include "grid_synchronizer.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NeighborhoodMaxCriterion::NeighborhoodMaxCriterion(
    Model & model, const ElementTypeMapReal & quad_coordinates,
    const ID & criterion_id, const ID & id)
    : NeighborhoodBase(model, quad_coordinates, id),
      Parsable(ParserType::_non_local, id), is_highest("is_highest", id),
      criterion(criterion_id, id) {

  AKANTU_DEBUG_IN();

  this->registerParam("radius", neighborhood_radius, 100.,
                      _pat_parsable | _pat_readable, "Non local radius");

  auto & mesh = this->model.getMesh();
  /// allocate the element type map arrays for _not_ghosts: One entry per quad
  GhostType ghost_type = _not_ghost;
  for (auto type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto new_size = this->quad_coordinates(type, ghost_type).size();
    this->is_highest.alloc(new_size, 1, type, ghost_type, true);
    this->criterion.alloc(new_size, 1, type, ghost_type, 1.);
  }

  /// criterion needs allocation also for ghost
  ghost_type = _ghost;
  for (auto type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto new_size = this->quad_coordinates(type, ghost_type).size();
    this->criterion.alloc(new_size, 1, type, ghost_type, true);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NeighborhoodMaxCriterion::~NeighborhoodMaxCriterion() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::initNeighborhood() {
  AKANTU_DEBUG_IN();

  /// parse the input parameter
  const auto & parser = getStaticParser();
  const auto & section_neighborhood =
      *(parser.getSubSections(ParserType::_neighborhood).first);
  this->parseSection(section_neighborhood);

  AKANTU_DEBUG_INFO("Creating the grid");
  this->createGrid();

  /// insert the non-ghost quads into the grid
  this->insertAllQuads(_not_ghost);

  /// store the number of current ghost elements for each type in the mesh
  ElementTypeMap<Int> nb_ghost_protected;
  auto & mesh = this->model.getMesh();
  for (auto type : mesh.elementTypes(spatial_dimension, _ghost)) {
    nb_ghost_protected(mesh.getNbElement(type, _ghost), type, _ghost);
  }

  /// create the grid synchronizer
  this->createGridSynchronizer();

  /// insert the ghost quads into the grid
  this->insertAllQuads(_ghost);

  /// create the pair lists
  this->updatePairList();

  /// remove the unneccessary ghosts
  this->cleanupExtraGhostElements(nb_ghost_protected);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::createGridSynchronizer() {
  this->is_creating_grid = true;
  std::set<SynchronizationTag> tags;
  tags.insert(SynchronizationTag::_nh_criterion);

  std::stringstream sstr;
  sstr << id << ":grid_synchronizer";
  this->grid_synchronizer = std::make_unique<GridSynchronizer>(
      this->model.getMesh(), *spatial_grid, *this, tags, sstr.str(), false);
  this->is_creating_grid = false;
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::insertAllQuads(GhostType ghost_type) {
  IntegrationPoint q;
  q.ghost_type = ghost_type;
  Mesh & mesh = this->model.getMesh();

  for (auto type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    Int nb_element = mesh.getNbElement(type, ghost_type);
    Int nb_quad =
        this->model.getFEEngine().getNbIntegrationPoints(type, ghost_type);

    const Array<Real> & quads = this->quad_coordinates(type, ghost_type);
    q.type = type;

    auto quad = quads.begin(spatial_dimension);

    for (Int e = 0; e < nb_element; ++e) {
      q.element = e;
      for (Int nq = 0; nq < nb_quad; ++nq) {
        q.num_point = nq;
        q.global_num = q.element * nb_quad + nq;
        spatial_grid->insert(q, *quad);
        ++quad;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::findMaxQuads(
    std::vector<IntegrationPoint> & max_quads) {
  AKANTU_DEBUG_IN();

  /// clear the element type maps
  this->is_highest.zero();
  this->criterion.zero();

  /// update the values of the criterion
  this->model.updateDataForNonLocalCriterion(criterion);

  /// start the exchange the value of the criterion on the ghost elements
  this->asynchronousSynchronize(SynchronizationTag::_nh_criterion);

  /// compare to not-ghost neighbors
  checkNeighbors(_not_ghost);

  /// finish the exchange
  this->waitEndSynchronize(SynchronizationTag::_nh_criterion);

  /// compare to ghost neighbors
  checkNeighbors(_ghost);

  /// extract the quads with highest criterion in their neighborhood
  IntegrationPoint quad;
  quad.ghost_type = _not_ghost;
  Mesh & mesh = this->model.getMesh();

  for (auto type : mesh.elementTypes(spatial_dimension, _not_ghost)) {
    quad.type = type;
    Int nb_quadrature_points =
        this->model.getFEEngine().getNbIntegrationPoints(type, _not_ghost);

    /// loop over is_highest for the current element type
    for (auto data : enumerate(is_highest(type, _not_ghost))) {
      const auto & is_highest = std::get<1>(data);
      if (is_highest) {
        auto q = std::get<0>(data);
        /// gauss point has the highest stress in his neighbourhood
        quad.element = q / nb_quadrature_points;
        quad.global_num = q;
        quad.num_point = q % nb_quadrature_points;
        max_quads.push_back(quad);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::checkNeighbors(GhostType ghost_type2) {
  AKANTU_DEBUG_IN();

  // Compute the weights
  for (auto & pair : pair_list[ghost_type2]) {
    const auto & lq1 = pair.first;
    const auto & lq2 = pair.second;

    Array<bool> & has_highest_eq_stress_1 =
        is_highest(lq1.type, lq1.ghost_type);

    const Array<Real> & criterion_1 = this->criterion(lq1.type, lq1.ghost_type);
    const Array<Real> & criterion_2 = this->criterion(lq2.type, lq2.ghost_type);

    if (criterion_1(lq1.global_num) < criterion_2(lq2.global_num)) {
      has_highest_eq_stress_1(lq1.global_num) = false;
    } else if (ghost_type2 != _ghost) {
      Array<bool> & has_highest_eq_stress_2 =
          is_highest(lq2.type, lq2.ghost_type);
      has_highest_eq_stress_2(lq2.global_num) = false;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NeighborhoodMaxCriterion::cleanupExtraGhostElements(
    const ElementTypeMap<Int> & nb_ghost_protected) {

  auto & mesh = this->model.getMesh();
  /// create remove elements event
  RemovedElementsEvent remove_elem(mesh);
  /// create set of ghosts to keep
  std::set<Element> relevant_ghost_elements;
  for (auto & pair : pair_list[_ghost]) {
    const auto & q2 = pair.second;
    relevant_ghost_elements.insert(q2);
  }

  Array<Element> ghosts_to_erase(0);

  Element element;
  element.ghost_type = _ghost;
  auto end = relevant_ghost_elements.end();
  for (const auto & type : mesh.elementTypes(spatial_dimension, _ghost)) {
    element.type = type;
    auto nb_ghost_elem = mesh.getNbElement(type, _ghost);
    decltype(nb_ghost_elem) nb_ghost_elem_protected = 0;
    try {
      nb_ghost_elem_protected = nb_ghost_protected(type, _ghost);
    } catch (...) {
    }

    if (!remove_elem.getNewNumbering().exists(type, _ghost)) {
      remove_elem.getNewNumbering().alloc(nb_ghost_elem, 1, type, _ghost);
    } else {
      remove_elem.getNewNumbering(type, _ghost).resize(nb_ghost_elem);
    }
    auto & new_numbering = remove_elem.getNewNumbering(type, _ghost);
    for (Int g = 0; g < nb_ghost_elem; ++g) {
      element.element = g;
      if (element.element >= nb_ghost_elem_protected &&
          relevant_ghost_elements.find(element) == end) {
        ghosts_to_erase.push_back(element);
        new_numbering(element.element) = -1;
      }
    }
    /// renumber remaining ghosts
    Int ng = 0;
    for (Int g = 0; g < nb_ghost_elem; ++g) {
      if (new_numbering(g) != Int(-1)) {
        new_numbering(g) = ng;
        ++ng;
      }
    }
  }

  mesh.sendEvent(remove_elem);
  this->onElementsRemoved(ghosts_to_erase, remove_elem.getNewNumbering(),
                          remove_elem);
}

} // namespace akantu
