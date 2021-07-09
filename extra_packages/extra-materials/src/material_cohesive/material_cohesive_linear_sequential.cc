/**
 * @file   material_cohesive_linear.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
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
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
/* ------------------------------------------------------------------ */
#include "material_cohesive_linear_sequential.hh"
/* ------------------------------------------------------------------ */

namespace akantu {

/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
MaterialCohesiveLinearSequential<spatial_dimension>::
    MaterialCohesiveLinearSequential(SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<spatial_dimension>(model, id),
      scalar_tractions("scalar_tractions", *this),
      normal_tractions("normal_tractions", *this),
      effective_stresses("effective_stresses", *this) {
  this->registerParam("delta_deviation", delta_deviation, 0., _pat_parsmod,
                      "Range within which quad points will be damaged");
}
/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  scalar_tractions.initialize(1);
  normal_tractions.initialize(spatial_dimension);
  effective_stresses.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::checkInsertion(
    bool check_only) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh_facets = this->model->getMeshFacets();

  for (auto && type_facet : mesh_facets.elementTypes(spatial_dimension - 1)) {

    // find the facet with the highest tensile stress
    auto max_stress_data = findCriticalFacet(type_facet);
    auto max_stress_facet = std::get<0>(max_stress_data);
    auto max_stress = std::get<1>(max_stress_data);
    auto crack_number = std::get<2>(max_stress_data);

    // communicate between processors the highest stress
    Real local_max_stress = max_stress;
    auto && comm = akantu::Communicator::getWorldCommunicator();
    comm.allReduce(max_stress, SynchronizerOperation::_max);

    // if the effective max stress is 0 or < 1 ->skip
    if (max_stress < 1)
      continue;

    // if the max stress at this proc != global max stress -> skip
    if (local_max_stress != max_stress)
      continue;

    // duplicate single facet and insert cohesive element
    std::map<UInt, UInt> facet_nb_crack_nb;
    facet_nb_crack_nb[max_stress_facet] = crack_number;
    insertCohesiveElements(facet_nb_crack_nb, type_facet, check_only);
  }

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(el_type);

  /// define iterators
  auto traction_it =
      this->tractions(el_type, ghost_type).begin(spatial_dimension);
  auto opening_it = this->opening(el_type, ghost_type).begin(spatial_dimension);
  auto contact_traction_it =
      this->contact_tractions(el_type, ghost_type).begin(spatial_dimension);
  auto contact_opening_it =
      this->contact_opening(el_type, ghost_type).begin(spatial_dimension);
  auto normal_opening_norm_it =
      this->normal_opening_norm(el_type, ghost_type).begin();
  auto tangential_opening_norm_it =
      this->tangential_opening_norm(el_type, ghost_type).begin();

  auto normal_it = normal.begin(spatial_dimension);
  auto traction_end =
      this->tractions(el_type, ghost_type).end(spatial_dimension);
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto insertion_stress_it =
      this->insertion_stress(el_type, ghost_type).begin(spatial_dimension);
  auto && element_filter = this->element_filter(el_type, ghost_type);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  // if (not this->model->isDefaultSolverExplicit())
  //   this->delta_max(el_type, ghost_type)
  //       .copy(this->delta_max.previous(el_type, ghost_type));

  /// loop on each quadrature point
  for (UInt i = 0; traction_it != traction_end; ++traction_it, ++opening_it,
            ++normal_it, ++sigma_c_it, ++delta_max_it, ++delta_c_it,
            ++damage_it, ++contact_traction_it, ++insertion_stress_it,
            ++contact_opening_it, ++normal_opening_norm_it,
            ++tangential_opening_norm_it, ++i) {
    bool penetration{false};
    UInt el_nb = floor(i / nb_quad_cohesive);
    __attribute__((unused)) auto critical_coh_el = element_filter(el_nb);
    this->computeTractionOnQuad(
        *traction_it, *opening_it, *normal_it, *delta_max_it, *delta_c_it,
        *insertion_stress_it, *sigma_c_it, normal_opening, tangential_opening,
        *normal_opening_norm_it, *tangential_opening_norm_it, *damage_it,
        penetration, *contact_traction_it, *contact_opening_it);
  }

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<
    spatial_dimension>::computeTangentTraction(ElementType el_type,
                                               Array<Real> & tangent_matrix,
                                               const Array<Real> & normal,
                                               GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(el_type);

  /// define iterators
  auto tangent_it = tangent_matrix.begin(spatial_dimension, spatial_dimension);

  auto tangent_end = tangent_matrix.end(spatial_dimension, spatial_dimension);

  auto normal_it = normal.begin(spatial_dimension);

  auto normal_opening_norm_it =
      this->normal_opening_norm(el_type, ghost_type).begin();
  auto tangential_opening_norm_it =
      this->tangential_opening_norm(el_type, ghost_type).begin();

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto sigma_c_it = this->sigma_c_eff(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto && element_filter = this->element_filter(el_type, ghost_type);

  for (UInt i = 0; tangent_it != tangent_end; ++tangent_it, ++normal_it,
            ++delta_max_it, ++sigma_c_it, ++delta_c_it, ++damage_it,
            ++normal_opening_norm_it, ++tangential_opening_norm_it, ++i) {
    UInt el_nb = floor(i / nb_quad_cohesive);
    __attribute__((unused)) auto critical_coh_el = element_filter(el_nb);

    this->computeTangentTractionOnQuad(
        *tangent_it, *delta_max_it, *delta_c_it, *sigma_c_it, *normal_it,
        *normal_opening_norm_it, *tangential_opening_norm_it, *damage_it);
  }

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------- */
template <UInt spatial_dimension>
UInt MaterialCohesiveLinearSequential<spatial_dimension>::updateDeltaMax(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(el_type);

  /// define iterators
  auto normal_opening_norm_it =
      this->normal_opening_norm(el_type, ghost_type).begin();
  auto normal_opening_norm_end =
      this->normal_opening_norm(el_type, ghost_type).end();
  auto tangential_opening_norm_it =
      this->tangential_opening_norm(el_type, ghost_type).begin();

  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto && element_filter = this->element_filter(el_type, ghost_type);

  Vector<Real> normal_opening(spatial_dimension);
  Vector<Real> tangential_opening(spatial_dimension);

  /// loop on each quadrature point
  UInt nb_updated_quads(0);
  for (UInt i = 0; normal_opening_norm_it != normal_opening_norm_end;
       ++normal_opening_norm_it, ++tangential_opening_norm_it, ++damage_it,
            ++delta_max_it, ++delta_c_it, ++i) {
    UInt el_nb = floor(i / nb_quad_cohesive);
    __attribute__((unused)) auto critical_coh_el = element_filter(el_nb);
    auto delta_max_increased = this->updateDeltaMaxOnQuad(
        *normal_opening_norm_it, *tangential_opening_norm_it, *damage_it,
        *delta_max_it, *delta_c_it);
    if (delta_max_increased)
      nb_updated_quads++;
  }

  AKANTU_DEBUG_OUT();
  return nb_updated_quads;
}

/* ------------------------------------------------------------------- */
template <UInt spatial_dimension>
std::tuple<Real, Element, UInt>
MaterialCohesiveLinearSequential<spatial_dimension>::computeMaxDeltaMaxExcess(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(el_type);

  /// define iterators
  auto normal_opening_norm_it =
      this->normal_opening_norm(el_type, ghost_type).begin();
  auto normal_opening_norm_end =
      this->normal_opening_norm(el_type, ghost_type).end();
  auto tangential_opening_norm_it =
      this->tangential_opening_norm(el_type, ghost_type).begin();

  auto delta_max_it = this->delta_max(el_type, ghost_type).begin();
  auto delta_c_it = this->delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = this->damage(el_type, ghost_type).begin();
  auto && element_filter = this->element_filter(el_type, ghost_type);

  Real max_delta_max_excess{0};
  /// loop on each quadrature point
  Element critical_coh_el{ElementNull};
  UInt nb_penetr_quads{0};
  for (UInt i = 0; normal_opening_norm_it != normal_opening_norm_end;
       ++normal_opening_norm_it, ++tangential_opening_norm_it, ++damage_it,
            ++delta_max_it, ++delta_c_it, ++i) {
    UInt el_nb = floor(i / nb_quad_cohesive);
    bool penetration{false};
    auto delta_max_excess = this->computeDeltaMaxExcessOnQuad(
        *normal_opening_norm_it, *tangential_opening_norm_it, *damage_it,
        *delta_max_it, *delta_c_it, penetration);
    if (delta_max_excess > max_delta_max_excess) {
      max_delta_max_excess = delta_max_excess;
      critical_coh_el.element = element_filter(el_nb);
      critical_coh_el.type = el_type;
      critical_coh_el.ghost_type = ghost_type;
    }
    if (penetration)
      nb_penetr_quads++;
  }

  AKANTU_DEBUG_OUT();
  return std::make_tuple(max_delta_max_excess, critical_coh_el,
                         nb_penetr_quads);
}
/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
std::tuple<UInt, Real, UInt>
MaterialCohesiveLinearSequential<spatial_dimension>::findCriticalFacet(
    const ElementType & type_facet) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  // MeshUtils::fillElementToSubElementsData(mesh_facets);
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, spatial_dimension).begin();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  Array<Element> candidate_facets;
  Array<UInt> coh_crack_nbs;
  Array<Real> highest_stresses;
  Array<UInt> coh_neighbors_nb;
  std::set<Element> crack_contour;
  auto output =
      std::make_tuple(UInt(-1), std::numeric_limits<Real>::min(), UInt(-1));

  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  const auto & facets_check = inserter.getCheckFacets(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  auto & scal_tractions = scalar_tractions(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & f_stress = this->model->getStressOnFacets(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & eff_stresses = effective_stresses(type_facet);

  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(type_cohesive);

  AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                      "The cohesive element and the corresponding facet do "
                      "not have the same numbers of integration points");
#endif
  // skip if no facets of this type are present
  UInt nb_facet = f_filter.size();
  if (nb_facet == 0)
    return output;

  Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
  UInt sp2 = spatial_dimension * spatial_dimension;

  const auto & tangents = this->model->getTangents(type_facet);
  const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                             .getNormalsOnIntegrationPoints(type_facet);
  auto normal_begin = normals.begin(spatial_dimension);
  auto tangent_begin = tangents.begin(tangents.getNbComponent());
  auto facet_stress_begin =
      f_stress.begin(spatial_dimension, spatial_dimension * 2);

  // loop over each facet belonging to this material
  for (auto && data :
       zip(f_filter, sigma_limits, eff_stresses,
           make_view(scal_tractions, nb_quad_facet),
           make_view(norm_tractions, spatial_dimension, nb_quad_facet))) {
    auto facet_nb = std::get<0>(data);
    Element facet{type_facet, facet_nb, _not_ghost};
    auto & sigma_limit = std::get<1>(data);
    auto & eff_stress = std::get<2>(data);
    auto & stress_check = std::get<3>(data);
    auto & normal_traction = std::get<4>(data);

    // skip facets where check shouldn't be inserted (or already inserted)
    if (!facets_check(facet_nb))
      continue;

    // compute the effective norm on each quadrature point of the facet
    for (UInt q : arange(nb_quad_facet)) {
      UInt current_quad = facet_nb * nb_quad_facet + q;
      const Vector<Real> & normal = normal_begin[current_quad];
      const Vector<Real> & tangent = tangent_begin[current_quad];
      const Matrix<Real> & facet_stress_it = facet_stress_begin[current_quad];

      // compute average stress on the current quadrature point
      Matrix<Real> stress_1(facet_stress_it.storage(), spatial_dimension,
                            spatial_dimension);

      Matrix<Real> stress_2(facet_stress_it.storage() + sp2, spatial_dimension,
                            spatial_dimension);

      stress_tmp.copy(stress_1);
      stress_tmp += stress_2;
      stress_tmp /= 2.;

      Vector<Real> normal_traction_vec(normal_traction(q));

      // compute normal and effective stress
      stress_check(q) = this->computeEffectiveNorm(stress_tmp, normal, tangent,
                                                   normal_traction_vec);
    }

    // verify if the effective stress overcomes the threshold
    Real final_stress = stress_check.mean();

    // normalize by the limit stress and skip non-stressed facets
    eff_stress = final_stress / sigma_limit;
    if (eff_stress < 1)
      continue;

    // check to how many cohesive elements this facet is connected
    const Vector<Element> & sub_to_facet =
        mesh_facets.getSubelementToElement(facet);
    std::vector<std::set<Element>> coh_neighbors(sub_to_facet.size());

    for (auto && data : enumerate(sub_to_facet)) {
      auto & subel = std::get<1>(data);
      auto & i = std::get<0>(data);
      auto & connected_elements_dim_min_1 =
          mesh_facets.getElementToSubelement(subel);
      for (auto & connected_element_dim_min_1 : connected_elements_dim_min_1) {
        auto & connected_elements_dim =
            mesh_facets.getElementToSubelement(connected_element_dim_min_1);
        for (auto & connected_element_dim : connected_elements_dim) {
          if (mesh.getKind(connected_element_dim.type) == _ek_cohesive) {
            coh_neighbors[i].emplace(connected_element_dim);
            crack_contour.insert(subel);
          }
        }
      }
    }

    // see if a single tip crack is present
    Array<UInt> single_tip_subs;
    for (UInt i : arange(sub_to_facet.size())) {
      if (coh_neighbors[i].size() == 1) {
        single_tip_subs.push_back(i);
        // experimental feature (all neighbors have to be single tips)
      } else if (coh_neighbors[i].size() > 1) {
        goto endloop;
      }
    }

    // if no coh els are connected - skip facet
    if (not single_tip_subs.size())
      goto endloop;

    // compute distances between barycenters and sharp angles between
    // facets; condition distance and angle
    UInt potential_crack_nb;
    for (auto tip_sub : single_tip_subs) {
      auto coh_el = *coh_neighbors[tip_sub].begin();
      auto & crack_numbers =
          mesh.getData<UInt>("crack_numbers", coh_el.type, coh_el.ghost_type);
      potential_crack_nb = crack_numbers(coh_el.element);

      // get inscribed diameter
      Real facet_indiam =
          MeshUtils::getInscribedCircleDiameter(*(this->model), facet);

      // get distance between two barycenters
      auto facet_to_coh_el = mesh_facets.getSubelementToElement(coh_el)(0);
      auto dist = MeshUtils::distanceBetweenIncentersCorrected(
          mesh_facets, facet, facet_to_coh_el);

      // ad-hoc rule on barycenters spacing
      // it should discard all elements under sharp angle
      if (dist < 0.8 * facet_indiam)
        goto endloop;

      // find facets bounding cohesive element
      const Vector<Element> & facets_to_cohesive =
          mesh_facets.getSubelementToElement(coh_el);
      // compute abs(dot) product between 2 normals & discard sharp angle
      Real dot = MeshUtils::cosSharpAngleBetween2Facets(*(this->model), facet,
                                                        facets_to_cohesive(0));

      if (dot < 0.7)
        goto endloop;
    }

    // add a candidate facet into a pool
    candidate_facets.push_back(facet);
    coh_crack_nbs.push_back(potential_crack_nb);
    highest_stresses.push_back(eff_stress);
    coh_neighbors_nb.push_back(single_tip_subs.size());

  endloop:;
  }
  if (not candidate_facets.size())
    return output;

  // insert element with the highest stress
  std::map<Real, UInt> stress_map;
  for (UInt i : arange(highest_stresses.size())) {
    stress_map.emplace(highest_stresses(i), i);
  }
  auto critical_facet_pos = stress_map.rbegin()->second;
  output = std::make_tuple(candidate_facets(critical_facet_pos).element,
                           highest_stresses(critical_facet_pos),
                           coh_crack_nbs(critical_facet_pos));

  return output;
}

/* ----------------------------------------------------------------- */

template <UInt spatial_dimension>
std::map<UInt, UInt> MaterialCohesiveLinearSequential<spatial_dimension>::
    findCriticalFacetsOnContour(
        const std::map<Element, Element> & contour_subfacets_coh_el,
        const std::map<Element, UInt> & surface_subfacets_crack_nb,
        const std::set<UInt> & contour_nodes,
        const std::set<UInt> & surface_nodes) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  UInt current_mat_index = this->model->getMaterialIndex(this->getName());
  // map is necessary to keep facet nbs sorted (otherwise insertion is messed
  // up)
  std::map<UInt, UInt> facet_nbs_crack_nbs;
  std::set<UInt> facet_nbs;
  std::set<Element> visited_facets;

  // TODO : make iteration on multiple facet types
  auto type_facet = *mesh_facets.elementTypes(spatial_dimension - 1).begin();
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  // skip if no facets of this type are present
  if (not f_filter.size()) {
    return facet_nbs_crack_nbs;
  }

  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);
#ifndef AKANTU_NDEBUG
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(type_cohesive);

  AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                      "The cohesive element and the corresponding facet do "
                      "not have the same numbers of integration points");
#endif
  CohesiveElementInserter & inserter = this->model->getElementInserter();
  const auto & facets_check = inserter.getCheckFacets(type_facet);
  auto & scal_tractions = scalar_tractions(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & f_stress = this->model->getStressOnFacets(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & eff_stresses = effective_stresses(type_facet);
  const auto & tangents = this->model->getTangents(type_facet);
  const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                             .getNormalsOnIntegrationPoints(type_facet);
  auto scal_tractions_it =
      scal_tractions.begin_reinterpret(nb_quad_facet, f_filter.size());
  auto norm_tractions_it = norm_tractions.begin_reinterpret(
      spatial_dimension, nb_quad_facet, f_filter.size());
  auto normal_it = normals.begin(spatial_dimension);
  auto tangent_it = tangents.begin(tangents.getNbComponent());
  auto facet_stress_it =
      f_stress.begin(spatial_dimension, spatial_dimension * 2);

  Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
  UInt sp2 = spatial_dimension * spatial_dimension;

  for (auto && pair : contour_subfacets_coh_el) {
    auto & contour_el = pair.first;
    auto & coh_el = pair.second;

    Real max_stress = std::numeric_limits<Real>::min();
    Element max_stress_facet{ElementNull};
    auto & facets_to_contour = mesh_facets.getElementToSubelement(contour_el);
    // identify the crack number and cohesive nodes
    UInt crack_nb = mesh.getData<UInt>("crack_numbers", coh_el.type,
                                       coh_el.ghost_type)(coh_el.element);

    auto && coh_facets_vec = mesh_facets.getSubelementToElement(coh_el);
    Array<Element> coh_facets;
    for (auto coh_facet : coh_facets_vec) {
      auto && subfacets = mesh_facets.getSubelementToElement(coh_facet);
      auto ret = std::find(subfacets.storage(),
                           subfacets.storage() + subfacets.size(), contour_el);
      if (ret == subfacets.storage() + subfacets.size())
        continue;
      coh_facets.push_back(coh_facet);
    }

    Vector<UInt> cohesive_nodes = mesh.getConnectivity(coh_el);
    AKANTU_DEBUG_ASSERT(coh_facets.size(),
                        "Cohesive element on countour not identified");

    // proceed to identifying the potential insertion facet
    for (auto & facet : facets_to_contour) {
      // check if the facet was already visited
      auto ret = visited_facets.emplace(facet);
      if (not ret.second)
        continue;

      // skip ghost facets
      if (facet.ghost_type == _ghost) {
        continue;
      }

      // skip cohesive facets
      if (coh_facets.find(facet) != UInt(-1)) {
        continue;
      }

      // check if the facet belongs to this material
      auto & facet_mat_index = this->model->getFacetMaterial(
          facet.type, facet.ghost_type)(facet.element);
      if (facet_mat_index != current_mat_index) {
        continue;
      }

      // check on check_facets
      if (not facets_check(facet.element)) {
      nextfacet:;
        continue;
      }

      // discard facets intersecting crack surface (unless different cracks)
      auto && subfacet_to_facet = mesh_facets.getSubelementToElement(facet);
      for (auto & subfacet : subfacet_to_facet) {
        auto ret = surface_subfacets_crack_nb.find(subfacet);
        if (ret != surface_subfacets_crack_nb.end()) {
          if (crack_nb == ret->second) {
            goto nextfacet;
          }
        }
      }

      // discard facets touching the surface
      Vector<UInt> facet_nodes = mesh_facets.getConnectivity(facet);
      for (auto & facet_node : facet_nodes) {
        if (surface_nodes.find(facet_node) != surface_nodes.end()) {
          goto nextfacet;
        }
      }

      // discard facets touching the contour
      UInt nb_subfacets_on_contour(0);
      for (auto & subfacet : subfacet_to_facet) {
        if (contour_subfacets_coh_el.find(subfacet) !=
            contour_subfacets_coh_el.end()) {
          nb_subfacets_on_contour++;
        }
      }
      if (nb_subfacets_on_contour == 1) {
        UInt nb_nodes_on_contour(0);
        Vector<UInt> facet_nodes = mesh_facets.getConnectivity(facet);
        for (auto & facet_node : facet_nodes) {
          if (contour_nodes.find(facet_node) != contour_nodes.end()) {
            nb_nodes_on_contour++;
          }
        }
        if (nb_nodes_on_contour == 3) {
          goto nextfacet;
        }
      }

      // discard facets under sharp angle with crack
      Real facet_indiam =
          MeshUtils::getInscribedCircleDiameter(*(this->model), facet);
      auto dist = MeshUtils::distanceBetweenIncentersCorrected(
          mesh_facets, facet, coh_facets[0]);
      // ad-hoc rule on barycenters spacing
      // it should discard all elements under sharp angle
      if (dist < 0.8 * facet_indiam) {
        continue;
      }

      // compute abs(dot) product between 2 normals & discard sharp angle
      Real dot = MeshUtils::cosSharpAngleBetween2Facets(*(this->model), facet,
                                                        coh_facets[0]);
      if (dot < 0.8) {
        continue;
      }

      // // non-local check on facet normal (all nearby facets)
      // CSR<Element> nodes_to_facets;
      // MeshUtils::buildNode2Elements(mesh, nodes_to_facets,
      //                               spatial_dimension - 1, _ek_regular);
      // Vector<Real> facet_normal = normal_it[facet.element * nb_quad_facet];
      // Vector<Real> final_normal(spatial_dimension);
      // std::set<Element> visited_facets;
      // Vector<Real> ref_vector(spatial_dimension);
      // for (auto node : cohesive_nodes) {
      //   for (auto & near_facet : nodes_to_facets.getRow(node)) {
      //     if (not visited_facets.size()) {
      //       const auto & near_facets_normals =
      //           this->model->getFEEngine("FacetsFEEngine")
      //               .getNormalsOnIntegrationPoints(near_facet.type,
      //                                              near_facet.ghost_type);
      //       auto near_facet_normal_it =
      //           near_facets_normals.begin(spatial_dimension);

      //       ref_vector =
      //           near_facet_normal_it[near_facet.element * nb_quad_facet];
      //     }
      //     auto ret = visited_facets.emplace(near_facet);
      //     if (ret.second) {
      //       const auto & near_facets_normals =
      //           this->model->getFEEngine("FacetsFEEngine")
      //               .getNormalsOnIntegrationPoints(near_facet.type,
      //                                              near_facet.ghost_type);
      //       auto near_facet_normal_it =
      //           near_facets_normals.begin(spatial_dimension);
      //       Vector<Real> near_facet_normal =
      //           near_facet_normal_it[near_facet.element * nb_quad_facet];
      //       Real test_dot = ref_vector.dot(near_facet_normal);
      //       if (test_dot < 0)
      //         near_facet_normal *= -1;

      //       final_normal += near_facet_normal;
      //     }
      //   }
      // }
      // final_normal /= final_normal.norm();

      // // compute abs(dot) product between 2 normals & discard sharp angle
      // Real dot = facet_normal.dot(final_normal);
      // dot = std::abs(dot);

      // if (dot < 0.8)
      //   continue;

      // check stress
      auto facet_local_nb = f_filter.find(facet.element);
      auto & sigma_limit = sigma_limits(facet_local_nb);
      auto & eff_stress = eff_stresses(facet_local_nb);
      Vector<Real> stress_check(scal_tractions_it[facet_local_nb]);
      Matrix<Real> normal_traction(norm_tractions_it[facet_local_nb]);

      // compute the effective norm on each quadrature point of the facet
      for (UInt q : arange(nb_quad_facet)) {
        UInt current_quad = facet.element * nb_quad_facet + q;
        const Vector<Real> & normal = normal_it[current_quad];
        const Vector<Real> & tangent = tangent_it[current_quad];
        const Matrix<Real> & facet_stress = facet_stress_it[current_quad];

        // compute average stress on the current quadrature point
        Matrix<Real> stress_1(facet_stress.storage(), spatial_dimension,
                              spatial_dimension);

        Matrix<Real> stress_2(facet_stress.storage() + sp2, spatial_dimension,
                              spatial_dimension);

        stress_tmp.copy(stress_1);
        stress_tmp += stress_2;
        stress_tmp /= 2.;

        Vector<Real> normal_traction_vec(normal_traction(q));

        // compute normal and effective stress
        stress_check(q) = this->computeEffectiveNorm(
            stress_tmp, normal, tangent, normal_traction_vec);
      }

      // verify if the effective stress overcomes the threshold
      Real final_stress = stress_check.mean();

      // normalize by the limit stress and skip non-stressed facets
      eff_stress = final_stress / sigma_limit;

      if (eff_stress < 1)
        continue;

      // determine the most stressed facet connected to this contour
      if (eff_stress > max_stress) {
        max_stress = eff_stress;
        max_stress_facet = facet;
      }
    }
    if (max_stress_facet != ElementNull) {
      auto ret = facet_nbs.emplace(max_stress_facet.element);
      if (ret.second) {
        facet_nbs_crack_nbs[max_stress_facet.element] = crack_nb;
      }
    }
  }
  return facet_nbs_crack_nbs;
}

/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<
    spatial_dimension>::computeEffectiveStresses() {

  AKANTU_DEBUG_IN();

  const Mesh & mesh_facets = this->model->getMeshFacets();

  for (auto && type_facet : mesh_facets.elementTypes(spatial_dimension - 1)) {

    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    auto & f_filter = this->facet_filter(type_facet);
    // skip if no facets of this type are present
    if (not f_filter.size()) {
      return;
    }

    UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                             .getNbIntegrationPoints(type_facet);
#ifndef AKANTU_NDEBUG
    UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif
    auto & scal_tractions = scalar_tractions(type_facet);
    auto & norm_tractions = normal_tractions(type_facet);
    const auto & f_stress = this->model->getStressOnFacets(type_facet);
    const auto & sigma_limits = this->sigma_c(type_facet);
    auto & eff_stresses = effective_stresses(type_facet);
    const auto & tangents = this->model->getTangents(type_facet);
    const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    auto scal_tractions_it =
        scal_tractions.begin_reinterpret(nb_quad_facet, f_filter.size());
    auto norm_tractions_it = norm_tractions.begin_reinterpret(
        spatial_dimension, nb_quad_facet, f_filter.size());
    auto normal_it = normals.begin(spatial_dimension);
    auto tangent_it = tangents.begin(tangents.getNbComponent());
    auto facet_stress_it =
        f_stress.begin(spatial_dimension, spatial_dimension * 2);

    Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
    UInt sp2 = spatial_dimension * spatial_dimension;

    for (auto && data : enumerate(f_filter)) {
      auto facet_local_nb = std::get<0>(data);
      auto facet_global_nb = std::get<1>(data);
      auto & sigma_limit = sigma_limits(facet_local_nb);
      auto & eff_stress = eff_stresses(facet_local_nb);
      Vector<Real> stress_check(scal_tractions_it[facet_local_nb]);
      Matrix<Real> normal_traction(norm_tractions_it[facet_local_nb]);

      // compute the effective norm on each quadrature point of the facet
      for (UInt q : arange(nb_quad_facet)) {
        UInt current_quad = facet_global_nb * nb_quad_facet + q;
        const Vector<Real> & normal = normal_it[current_quad];
        const Vector<Real> & tangent = tangent_it[current_quad];
        const Matrix<Real> & facet_stress = facet_stress_it[current_quad];

        // compute average stress on the current quadrature point
        Matrix<Real> stress_1(facet_stress.storage(), spatial_dimension,
                              spatial_dimension);

        Matrix<Real> stress_2(facet_stress.storage() + sp2, spatial_dimension,
                              spatial_dimension);

        stress_tmp.copy(stress_1);
        stress_tmp += stress_2;
        stress_tmp /= 2.;

        Vector<Real> normal_traction_vec(normal_traction(q));

        // compute normal and effective stress
        stress_check(q) = this->computeEffectiveNorm(
            stress_tmp, normal, tangent, normal_traction_vec);
      }

      // verify if the effective stress overcomes the threshold
      Real final_stress = stress_check.mean();

      // normalize by the limit stress and skip non-stressed facets
      eff_stress = final_stress / sigma_limit;
    }
  }
}

/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
std::map<UInt, UInt>
MaterialCohesiveLinearSequential<spatial_dimension>::findHolesOnContour(
    std::map<Element, Element> & contour_subfacets_coh_el,
    const std::map<Element, UInt> & /*surface_subfacets_crack_nb*/,
    const std::set<UInt> & /*surface_nodes*/) {
  AKANTU_DEBUG_IN();
  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  CohesiveElementInserter & inserter = this->model->getElementInserter();
  UInt current_mat_index = this->model->getMaterialIndex(this->getName());

  std::map<UInt, UInt> facet_nbs_crack_nbs;
  std::map<UInt, std::set<Element>> facet_contour_subfacets;
  std::map<UInt, UInt> facet_crack_nb;
  std::set<Element> facet_pool;

  // TODO : make iteration on multiple facet types
  auto type_facet = *mesh_facets.elementTypes(spatial_dimension - 1).begin();
  auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);
  const auto & facets_check = inserter.getCheckFacets(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  // auto & scal_tractions = scalar_tractions(type_facet);
  // auto & norm_tractions = normal_tractions(type_facet);
  // const auto & f_stress = this->model->getStressOnFacets(type_facet);
  // const auto & sigma_limits = this->sigma_c(type_facet);
  // auto & eff_stresses = effective_stresses(type_facet);
  // const auto & tangents = this->model->getTangents(type_facet);
  // const auto & normals = this->model->getFEEngine("FacetsFEEngine")
  //                            .getNormalsOnIntegrationPoints(type_facet);
  // auto normal_it = normals.begin(spatial_dimension);
  // auto scal_tractions_it =
  //     scal_tractions.begin_reinterpret(nb_quad_facet, f_filter.size());
  // auto norm_tractions_it = norm_tractions.begin_reinterpret(
  //     spatial_dimension, nb_quad_facet, f_filter.size());
  // auto tangent_it = tangents.begin(tangents.getNbComponent());
  // auto facet_stress_it =
  //     f_stress.begin(spatial_dimension, spatial_dimension * 2);

  // Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
  // UInt sp2 = spatial_dimension * spatial_dimension;

#ifndef AKANTU_NDEBUG
  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(type_cohesive);

  AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                      "The cohesive element and the corresponding facet do "
                      "not have the same numbers of integration points");
#endif
  // skip if no facets of this type are present
  UInt nb_facet = f_filter.size();
  if (nb_facet == 0)
    return facet_nbs_crack_nbs;

  for (auto && pair : contour_subfacets_coh_el) {
    auto & contour_el = pair.first;
    auto & coh_el = pair.second;
    auto & facets_to_contour = mesh_facets.getElementToSubelement(contour_el);
    // identify the crack number and cohesive nodes
    UInt crack_nb = mesh.getData<UInt>("crack_numbers", coh_el.type,
                                       coh_el.ghost_type)(coh_el.element);
    auto && coh_facets_vec = mesh_facets.getSubelementToElement(coh_el);
    Array<Element> coh_facets;
    for (auto coh_facet : coh_facets_vec) {
      auto && subfacets = mesh_facets.getSubelementToElement(coh_facet);
      auto ret = std::find(subfacets.storage(),
                           subfacets.storage() + subfacets.size(), contour_el);
      if (ret == subfacets.storage() + subfacets.size())
        continue;
      coh_facets.push_back(coh_facet);
    }
    AKANTU_DEBUG_ASSERT(coh_facets.size(),
                        "Cohesive element on countour not identified");

    // proceed to identifying the potential insertion facet
    for (auto & facet : facets_to_contour) {

      // skip ghost facets
      if (facet.ghost_type == _ghost) {
        continue;
      }

      // skip cohesive facets
      if (coh_facets.find(facet) != UInt(-1)) {
        continue;
      }

      // check if the facet belongs to this material
      auto & facet_mat_index = this->model->getFacetMaterial(
          facet.type, facet.ghost_type)(facet.element);
      if (facet_mat_index != current_mat_index) {
        continue;
      }
      // check check_facets
      if (not facets_check(facet.element)) {
        continue;
      }

      // // discard facets intersecting crack surface
      // auto && subfacet_to_facet = mesh_facets.getSubelementToElement(facet);
      // for (auto & subfacet : subfacet_to_facet) {
      //   if (surface_subfacets.find(subfacet) != surface_subfacets.end()) {
      //     goto nextfacet;
      //   }
      // }

      // // discard facets touching the surface
      // Vector<UInt> facet_nodes = mesh_facets.getConnectivity(facet);
      // for (auto & facet_node : facet_nodes) {
      //   if (surface_nodes.find(facet_node) != surface_nodes.end()) {
      //     goto nextfacet;
      //   }
      // }

      // discard facets under sharp angle with crack
      Real facet_indiam =
          MeshUtils::getInscribedCircleDiameter(*(this->model), facet);
      auto dist = MeshUtils::distanceBetweenIncentersCorrected(
          mesh_facets, facet, coh_facets[0]);
      // ad-hoc rule on barycenters spacing
      // it should discard all elements under sharp angle
      if (dist < 0.8 * facet_indiam) {
        continue;
      }

      // avoid facets under sharp angles with their cohesives
      Real dot = MeshUtils::cosSharpAngleBetween2Facets(*(this->model), facet,
                                                        coh_facets[0]);
      if (dot < 0.8)
        continue;

      // // discard not stressed or compressed facets
      // auto facet_local_nb = f_filter.find(facet.element);
      // auto & sigma_limit = sigma_limits(facet_local_nb);
      // auto & eff_stress = eff_stresses(facet_local_nb);
      // Vector<Real> stress_check(scal_tractions_it[facet_local_nb]);
      // Matrix<Real> normal_traction(norm_tractions_it[facet_local_nb]);

      // // compute the effective norm on each quadrature point of the facet
      // for (UInt q : arange(nb_quad_facet)) {
      //   UInt current_quad = facet.element * nb_quad_facet + q;
      //   const Vector<Real> & normal = normal_it[current_quad];
      //   const Vector<Real> & tangent = tangent_it[current_quad];
      //   const Matrix<Real> & facet_stress = facet_stress_it[current_quad];
      //   // compute average stress on the current quadrature point
      //   Matrix<Real> stress_1(facet_stress.storage(), spatial_dimension,
      //                         spatial_dimension);

      //   Matrix<Real> stress_2(facet_stress.storage() + sp2,
      //   spatial_dimension,
      //                         spatial_dimension);
      //   stress_tmp.copy(stress_1);
      //   stress_tmp += stress_2;
      //   stress_tmp /= 2.;

      //   Vector<Real> normal_traction_vec(normal_traction(q));
      //   // compute normal and effective stress
      //   stress_check(q) = this->computeEffectiveNorm(
      //       stress_tmp, normal, tangent, normal_traction_vec);
      // }

      // // verify if the effective stress overcomes the threshold
      // Real final_stress = stress_check.mean();

      // // normalize by the limit stress and skip non-stressed facets
      // eff_stress = final_stress / sigma_limit;
      // if (eff_stress <= 0)
      //   continue;

      // add facet into the map
      facet_contour_subfacets[facet.element].emplace(contour_el);
      facet_crack_nb[facet.element] = crack_nb;
      // add all passed facets into the pool
      facet_pool.emplace(facet);
    }
  }

  // get all facets having more than 1 connection to the contour
  // and clean facet_pool from the facets connected to these contours
  for (auto it = facet_contour_subfacets.begin();
       it != facet_contour_subfacets.end(); it++) {
    if (it->second.size() > 1) {
      facet_nbs_crack_nbs[it->first] = facet_crack_nb[it->first];

      // delete facets connected to above subfacets from the pool
      Element facet{type_facet, it->first, _not_ghost};
      for (auto & contour_subf : it->second) {
        auto & facets_to_taken_contour =
            mesh_facets.getElementToSubelement(contour_subf);
        for (auto & facet_to_taken_contour : facets_to_taken_contour) {
          if (facet_to_taken_contour == facet)
            continue;
          facet_pool.erase(facet_to_taken_contour);
        }
      }
    }
  }

  // determine all interconnected holes
  // initialize a graph
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
      Graph;
  Graph graph;
  // fill the graph in
  auto first = facet_pool.begin();
  for (auto it1 = first; it1 != facet_pool.end(); it1++) {
    // insert a single edge
    boost::add_edge(std::distance(first, it1), std::distance(first, it1),
                    graph);
    for (auto it2 = std::next(it1); it2 != facet_pool.end(); it2++) {
      auto data = MeshUtils::areFacetsConnected(mesh_facets, *it1, *it2);
      // facets are not connected
      if (not std::get<0>(data))
        continue;
      // connected but through the current crack tip
      for (auto & shared_sub : std::get<1>(data)) {
        if (contour_subfacets_coh_el.find(shared_sub) !=
            contour_subfacets_coh_el.end()) {
          goto endloop;
        }
      }
      // otherwise properly connected;
      boost::add_edge(std::distance(first, it1), std::distance(first, it2),
                      graph);
    endloop:;
    }
  }
  // connectivity and number of components of a graph
  std::vector<int> components(boost::num_vertices(graph));
  boost::connected_components(graph, &components[0]);
  if (not components.size()) {
    return facet_nbs_crack_nbs;
  }
  // compute component sizes
  std::map<int, UInt> components_size;
  for (auto component : components) {
    components_size[component]++;
  }
  // insert all components with size > 1
  for (auto it = first; it != facet_pool.end(); it++) {
    auto component = components[std::distance(first, it)];
    if (components_size[component] <= 1)
      continue;
    // otherwise insert this facet as it is a godly act
    facet_nbs_crack_nbs[it->element] = facet_crack_nb[it->element];
  }
  AKANTU_DEBUG_OUT();
  return facet_nbs_crack_nbs;
}

/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
std::tuple<std::map<Element, Element>, std::map<Element, UInt>, std::set<UInt>,
           std::set<UInt>>
MaterialCohesiveLinearSequential<spatial_dimension>::determineCrackSurface() {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  std::map<Element, Element> contour_subfacets_coh_el;
  std::map<Element, UInt> surface_subfacets_crack_nb;
  std::set<UInt> surface_nodes;
  std::set<UInt> contour_nodes;
  std::set<Element> visited_subfacets;

  // first loop on facets (doesn't include duplicated ones)
  for (auto & type_facet :
       mesh.elementTypes(spatial_dimension - 1, _not_ghost, _ek_regular)) {
    auto & filter = this->facet_filter(type_facet);
    if (filter.empty()) {
      continue;
    }

    for (auto & facet_nb : filter) {
      Element facet{type_facet, facet_nb, _not_ghost};
      auto && subfacets_to_facet = mesh_facets.getSubelementToElement(facet);
      for (auto & subfacet : subfacets_to_facet) {
        auto ret = visited_subfacets.emplace(subfacet);
        if (not ret.second)
          continue;

        std::set<Element> connected_cohesives;
        auto && facets_to_subfacet =
            mesh_facets.getElementToSubelement(subfacet);
        for (auto & facet : facets_to_subfacet) {
          auto & connected_to_facet = mesh_facets.getElementToSubelement(facet);
          for (auto & connected_el : connected_to_facet) {
            if (connected_el.type == _not_defined)
              goto nextfacet;
            if (mesh.getKind(connected_el.type) == _ek_cohesive) {
              connected_cohesives.emplace(connected_el);
            }
          }
        nextfacet:;
        }
        if (connected_cohesives.size() == 1) {
          contour_subfacets_coh_el[subfacet] = *connected_cohesives.begin();
          Vector<UInt> contour_subfacet_nodes =
              mesh_facets.getConnectivity(subfacet);
          for (auto & contour_subfacet_node : contour_subfacet_nodes) {
            contour_nodes.emplace(contour_subfacet_node);
          }
        } else if (connected_cohesives.size() > 1) {
          auto coh_element = *connected_cohesives.begin();
          UInt crack_nb =
              mesh.getData<UInt>("crack_numbers", coh_element.type,
                                 coh_element.ghost_type)(coh_element.element);
          surface_subfacets_crack_nb[subfacet] = crack_nb;
          // add subfacet nodes to the set
          for (auto & subfacet_node : mesh_facets.getConnectivity(subfacet))
            surface_nodes.emplace(subfacet_node);
        }
      }
    }
  }

  // second loop on cohesives (includes duplicated facets)
  for (auto & type_cohesive :
       mesh.elementTypes(spatial_dimension, _not_ghost, _ek_cohesive)) {
    auto & filter = this->element_filter(type_cohesive, _not_ghost);
    if (filter.empty()) {
      continue;
    }

    for (auto & cohesive_nb : filter) {
      Element cohesive{type_cohesive, cohesive_nb, _not_ghost};
      auto && facets_to_cohesive = mesh_facets.getSubelementToElement(cohesive);
      for (auto coh_facet : facets_to_cohesive) {
        auto && subfacets_to_facet =
            mesh_facets.getSubelementToElement(coh_facet);
        for (auto & subfacet : subfacets_to_facet) {
          auto ret = visited_subfacets.emplace(subfacet);
          if (not ret.second)
            continue;

          std::set<Element> connected_cohesives;
          auto && facets_to_subfacet =
              mesh_facets.getElementToSubelement(subfacet);
          for (auto & facet : facets_to_subfacet) {
            auto & connected_to_facet =
                mesh_facets.getElementToSubelement(facet);
            for (auto & connected_el : connected_to_facet) {
              if (connected_el.type == _not_defined)
                goto nextfacet2;
              if (mesh.getKind(connected_el.type) == _ek_cohesive) {
                connected_cohesives.emplace(connected_el);
              }
            }
          nextfacet2:;
          }
          if (connected_cohesives.size() == 1) {
            contour_subfacets_coh_el[subfacet] = *connected_cohesives.begin();
            Vector<UInt> contour_subfacet_nodes =
                mesh_facets.getConnectivity(subfacet);
            for (auto & contour_subfacet_node : contour_subfacet_nodes) {
              contour_nodes.emplace(contour_subfacet_node);
            }
          } else if (connected_cohesives.size() > 1) {
            auto coh_element = *connected_cohesives.begin();
            UInt crack_nb =
                mesh.getData<UInt>("crack_numbers", coh_element.type,
                                   coh_element.ghost_type)(coh_element.element);
            surface_subfacets_crack_nb[subfacet] = crack_nb;
            // add subfacet nodes to the set
            for (auto & subfacet_node : mesh_facets.getConnectivity(subfacet))
              surface_nodes.emplace(subfacet_node);
          }
        }
      }
    }
  }

  // remove external nodes from the surface nodes list
  for (auto & contour_node : contour_nodes) {
    surface_nodes.erase(contour_node);
  }

  return std::make_tuple(contour_subfacets_coh_el, surface_subfacets_crack_nb,
                         surface_nodes, contour_nodes);
}
/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::
    insertCohesiveElements(std::map<UInt, UInt> & facet_nbs_crack_nbs,
                           ElementType type_facet, bool check_only) {
  if (not facet_nbs_crack_nbs.size()) {
    return;
  }

  CohesiveElementInserter & inserter = this->model->getElementInserter();
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  auto & f_insertion = inserter.getInsertionFacets(type_facet);
  auto & scal_tractions = scalar_tractions(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & sig_c_eff = this->sigma_c_eff(type_cohesive);
  auto & del_c = this->delta_c_eff(type_cohesive);
  auto & del_max = this->delta_max(type_cohesive);
  auto & ins_stress = this->insertion_stress(type_cohesive);
  auto & trac_old = this->tractions.previous(type_cohesive);
  const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                             .getNormalsOnIntegrationPoints(type_facet);
  auto normal_begin = normals.begin(spatial_dimension);

  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);
  std::vector<Real> new_sigmas;
  std::vector<Vector<Real>> new_normal_traction;
  std::vector<Real> new_delta_c;
  std::vector<Real> new_delta_max;
  UInt nb_understressed{0};

  for (auto && pair : facet_nbs_crack_nbs) {
    auto & facet_nb = pair.first;
    // get facet's local id
    auto local_id = f_filter.find(facet_nb);
    AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                        "mismatch between global and local facet numbering"
                            << facet_nb);

    // mark the insertion of the cohesive element
    f_insertion(facet_nb) = true;

    if (check_only)
      continue;

    auto & sigma_limit = sigma_limits(local_id);
    Vector<Real> stress_check =
        make_view(scal_tractions, nb_quad_facet).begin()[local_id];
    Matrix<Real> normal_traction =
        make_view(norm_tractions, spatial_dimension, nb_quad_facet)
            .begin()[local_id];

    // store the new cohesive material parameters for each quad point
    for (UInt q = 0; q < nb_quad_facet; ++q) {
      UInt current_quad = facet_nb * nb_quad_facet + q;
      Real new_sigma = stress_check(q);
      Real new_del_max{0};
      Vector<Real> normal_traction_vec(normal_traction(q));
      const Vector<Real> & normal = normal_begin[current_quad];

      // workaround to insert understressed cohesives
      if (new_sigma < sigma_limit) {
        new_sigma = sigma_limit;
        normal_traction_vec = normal * sigma_limit;
        nb_understressed++;
      }

      // EXPERIMENTAL !!!!! FOR ESHELBY PROBLEM
      if (new_sigma > sigma_limit) {
        new_sigma = sigma_limit;
        normal_traction_vec = normal * sigma_limit;
      }

      Real new_delta;
      // set delta_c in function of G_c or a given delta_c value
      if (Math::are_float_equal(this->delta_c, 0.))
        new_delta = 2 * this->G_c / new_sigma;
      else
        new_delta = sigma_limit / new_sigma * this->delta_c;

      // artifitial way to give initial stiffness to cohesive (SLA)
      new_del_max = new_delta / 100;

      if (spatial_dimension != 3)
        normal_traction_vec *= -1.;

      new_sigmas.push_back(new_sigma);
      new_normal_traction.push_back(normal_traction_vec);
      new_delta_max.push_back(new_del_max);
      new_delta_c.push_back(new_delta);
    }
  }

  // update material data for the new elements
  UInt old_nb_quad_points = sig_c_eff.size();
  UInt new_nb_quad_points = new_sigmas.size();
  sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
  ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
  trac_old.resize(old_nb_quad_points + new_nb_quad_points);
  del_c.resize(old_nb_quad_points + new_nb_quad_points);
  del_max.resize(old_nb_quad_points + new_nb_quad_points);

  for (UInt q = 0; q < new_nb_quad_points; ++q) {
    sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
    del_c(old_nb_quad_points + q) = new_delta_c[q];
    del_max(old_nb_quad_points + q) = new_delta_max[q];
    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      ins_stress(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
      trac_old(old_nb_quad_points + q, dim) = new_normal_traction[q](dim);
    }
  }
  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(cohesive_linear_sequential,
                     MaterialCohesiveLinearSequential);

} // namespace akantu
