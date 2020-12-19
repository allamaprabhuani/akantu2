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
      normal_stresses("normal_stresses", *this),
      normal_tractions("normal_tractions", *this),
      effective_stresses("effective_stresses", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam(
      "insertion_threshold", insertion_threshold, Real(0.0),
      _pat_parsable | _pat_readable,
      "Portion of elements below the highest stress to be inserted");

  AKANTU_DEBUG_OUT();
}
/* ------------------------------------------------------------------ */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  normal_stresses.initialize(1);
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
    Array<std::pair<UInt, UInt>> facet_nb_crack_nb;
    facet_nb_crack_nb.push_back(std::make_pair(max_stress_facet, crack_number));
    insertCohesiveElements(facet_nb_crack_nb, type_facet, check_only);
  }

  AKANTU_DEBUG_OUT();
}
/* ----------------------------------------------------------------- */
template <>
void MaterialCohesiveLinearSequential<3>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();
  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, spatial_dimension).begin();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  Array<std::pair<UInt, UInt>> facet_nbs_crack_nbs;
  for (auto && facet_type : mesh_facets.elementTypes(2)) {

    Array<Element> candidate_facets;
    Array<UInt> coh_crack_nbs;
    Array<Real> highest_stresses;
    Array<UInt> coh_neighbors_nb;
    std::set<Element> crack_contour;

    ElementType type_cohesive = FEEngine::getCohesiveElementType(facet_type);
    const auto & facets_check = inserter.getCheckFacets(facet_type);
    auto & f_filter = this->facet_filter(facet_type);
    auto & norm_stresses = normal_stresses(facet_type);
    auto & norm_tractions = normal_tractions(facet_type);
    const auto & f_stress = this->model->getStressOnFacets(facet_type);
    const auto & sigma_limits = this->sigma_c(facet_type);
    auto & eff_stresses = effective_stresses(facet_type);

    UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                             .getNbIntegrationPoints(facet_type);

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
      continue;

    Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
    UInt sp2 = spatial_dimension * spatial_dimension;

    const auto & tangents = this->model->getTangents(facet_type);
    const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(facet_type);
    auto normal_begin = normals.begin(spatial_dimension);
    auto tangent_begin = tangents.begin(tangents.getNbComponent());
    auto facet_stress_begin =
        f_stress.begin(spatial_dimension, spatial_dimension * 2);

    // loop over each facet belonging to this material
    for (auto && data :
         zip(f_filter, sigma_limits, eff_stresses,
             make_view(norm_stresses, nb_quad_facet),
             make_view(norm_tractions, spatial_dimension, nb_quad_facet))) {
      auto facet_nb = std::get<0>(data);
      Element facet{facet_type, facet_nb, _not_ghost};
      auto & sigma_limit = std::get<1>(data);
      auto & eff_stress = std::get<2>(data);
      auto & stress_check = std::get<3>(data);
      auto & normal_traction = std::get<4>(data);

      // skip facets which shouldn't be inserted (or already inserted)
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

        Matrix<Real> stress_2(facet_stress_it.storage() + sp2,
                              spatial_dimension, spatial_dimension);

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

      // check to how many cohesive elements this facet is connected
      const Vector<Element> & sub_to_facet =
          mesh_facets.getSubelementToElement(facet);
      std::vector<std::set<Element>> coh_neighbors(sub_to_facet.size());

      for (auto && data : enumerate(sub_to_facet)) {
        auto & subel = std::get<1>(data);
        auto & i = std::get<0>(data);
        auto & connected_elements_dim_min_1 =
            mesh_facets.getElementToSubelement(subel);
        for (auto & connected_element_dim_min_1 :
             connected_elements_dim_min_1) {
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
      bool single_tip{false};
      Array<UInt> single_tip_subs;
      for (UInt i : arange(sub_to_facet.size())) {
        if (coh_neighbors[i].size() == 1) {
          single_tip = true;
          single_tip_subs.push_back(i);
          // experimental feature (all neighbors have to be single tips)
        } else if (coh_neighbors[i].size() > 1) {
          single_tip = false;
          break;
        }
      }

      // if no coh els are connected or no single tip present - skip face
      if (not single_tip)
        continue;

      // compute distances between barycenters and sharp angles between
      // facets; condition distance and angle
      bool sharp_angle{false};
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
        auto dist = MeshUtils::distanceBetween2Barycenters(mesh_facets, mesh,
                                                           facet, coh_el);

        // ad-hoc rule on barycenters spacing
        // it should discard all elements under sharp angle
        if (dist < 0.7 * facet_indiam) {
          sharp_angle = true;
          continue;
        }

        // find facets bounding cohesive element
        const Vector<Element> & facets_to_cohesive =
            mesh_facets.getSubelementToElement(coh_el);
        // compute abs(dot) product between 2 normals & discard sharp angle
        Real dot = MeshUtils::cosSharpAngleBetween2Facets(
            *(this->model), facet, facets_to_cohesive(0));

        if (dot < 0.7)
          sharp_angle = true;
      }

      if (sharp_angle)
        continue;

      // add a candidate facet into a pool
      candidate_facets.push_back(facet);
      coh_crack_nbs.push_back(potential_crack_nb);
      highest_stresses.push_back(eff_stress);
      coh_neighbors_nb.push_back(single_tip_subs.size());
    }

    if (not candidate_facets.size())
      continue;

    /*    // find the longest row of interconnected candidates
        // initialize a graph
        typedef boost::adjacency_list<boost::vecS, boost::vecS,
       boost::undirectedS> Graph; Graph graph;

        for (UInt i : arange(candidate_facets.size())) {
          // include single facet into graph by default
          boost::add_edge(i, i, graph);
          for (UInt j : arange(i + 1, candidate_facets.size())) {
            auto data = MeshUtils::areFacetsConnected(
                mesh_facets, candidate_facets(i), candidate_facets(j));
            // facets are not connected
            if (not std::get<0>(data))
              continue;
            // connected but through the current crack tip
            for (auto & shared_sub : std::get<1>(data)) {
              if (crack_contour.find(shared_sub) != crack_contour.end()) {
                goto endloop;
              }
            }
            // otherwise properly connected;
            boost::add_edge(i, j, graph);
          endloop:;
          }
        }

        // connectivity and number of components of a graph
        std::vector<int> component(boost::num_vertices(graph));
        boost::connected_components(graph, &component[0]);
        if (not component.size())
          continue;

        // create stress - facet map
        std::map<Real, int> stress_facet;
        for (int i : arange(highest_stresses.size())) {
          stress_facet[highest_stresses[i]] = i;
        }

        // identify highest stressed facet -> component
        int critical_component = component[stress_facet.rbegin()->second];

        // insert the longest row of connected facets
        for (UInt i : arange(component.size())) {
          if (component[i] == critical_component)
            facet_nbs_crack_nbs.push_back(
                std::make_pair(candidate_facets(i).element, coh_crack_nbs(i)));
        } */

    /*    // chose the most stressed element with highest nb of coh neighbors
        std::multimap<UInt, UInt> coh_neighbors_map;
        for (UInt i : arange(candidate_facets.size())) {
          coh_neighbors_map.emplace(coh_neighbors_nb(i), i);
        }
        // get highest nb of coh neighbors
        auto max_neighbors = coh_neighbors_map.rbegin()->first;
        // get range of facets with same nb of coh neighbors
        std::pair<std::multimap<UInt, UInt>::iterator,
                  std::multimap<UInt, UInt>::iterator>
            ret;
        ret = coh_neighbors_map.equal_range(max_neighbors);
        // chose the element with the highest stress
        std::map<Real, UInt> stress_map;
        for (std::multimap<UInt, UInt>::iterator it = ret.first; it !=
       ret.second;
             ++it) {
          stress_map.emplace(highest_stresses(it->second), it->second);
        }
        auto critical_facet_order = stress_map.rbegin()->second;
        facet_nbs_crack_nbs.push_back(
            std::make_pair(candidate_facets(critical_facet_order).element,
                           coh_crack_nbs(critical_facet_order)));
    */

    // insert element with the highest stress
    std::map<Real, UInt> stress_map;
    for (UInt i : arange(highest_stresses.size())) {
      stress_map.emplace(highest_stresses(i), i);
    }
    auto critical_facet_pos = stress_map.rbegin()->second;
    facet_nbs_crack_nbs.push_back(
        std::make_pair(candidate_facets(critical_facet_pos).element,
                       coh_crack_nbs(critical_facet_pos)));

    // insert cohesives
    insertCohesiveElements(facet_nbs_crack_nbs, facet_type, check_only);
  }
  AKANTU_DEBUG_OUT();
}

/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
bool MaterialCohesiveLinearSequential<spatial_dimension>::fillInHoles(
    ElementType facet_type, bool check_only) {
  AKANTU_DEBUG_IN();
  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  Array<std::pair<UInt, UInt>> facet_nbs_crack_nbs;
  Array<Element> understressed_candidates;
  Array<UInt> coh_crack_nbs;
  std::set<Element> crack_contour;

  ElementType type_cohesive = FEEngine::getCohesiveElementType(facet_type);
  const auto & facets_check = inserter.getCheckFacets(facet_type);
  auto & f_filter = this->facet_filter(facet_type);

  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(facet_type);

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
    return 0;

  // loop over each facet belonging to this material
  for (auto && facet_nb : f_filter) {
    Element facet{facet_type, facet_nb, _not_ghost};

    // skip facets which shouldn't be inserted (or already inserted)
    if (!facets_check(facet_nb))
      continue;

    Real facet_indiam =
        MeshUtils::getInscribedCircleDiameter(*(this->model), facet);

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
    UInt potential_crack_nb;
    std::map<UInt, UInt> crack_nb_repetitions;
    for (UInt i : arange(sub_to_facet.size())) {
      if (coh_neighbors[i].size() == 1) {
        // discard facets that form sharp angles with cohesives
        auto dist = MeshUtils::distanceBetween2Barycenters(
            mesh_facets, mesh, facet, *coh_neighbors[i].begin());
        if (dist < 0.7 * facet_indiam) {
          goto endloop;
        }
        single_tip_subs.push_back(i);
        // experimental feature (all neighbors have to be single tips)
      } else if (coh_neighbors[i].size() > 1) {
        goto endloop;
      }
    }
    if (single_tip_subs.size() < 2)
      goto endloop;

    // extract crack number
    for (auto tip_sub : single_tip_subs) {
      auto coh_el = *coh_neighbors[tip_sub].begin();
      auto & crack_numbers =
          mesh.getData<UInt>("crack_numbers", coh_el.type, coh_el.ghost_type);
      potential_crack_nb = crack_numbers(coh_el.element);
      ++crack_nb_repetitions[crack_numbers(coh_el.element)];
    }
    // see if any crack nb is met more than once
    for (auto it = crack_nb_repetitions.begin();
         it != crack_nb_repetitions.end(); it++) {
      if (it->second > 1)
        potential_crack_nb = it->first;
    }

    // fill in understressed candidates: low stress but many coh neighb
    understressed_candidates.push_back(facet);
    coh_crack_nbs.push_back(potential_crack_nb);

  endloop:;
  }

  if (not understressed_candidates.size())
    return false;

  // insert all understressed elements with multiple coh neighbors
  for (auto && data : zip(understressed_candidates, coh_crack_nbs)) {
    facet_nbs_crack_nbs.push_back(
        std::make_pair(std::get<0>(data).element, std::get<1>(data)));
  }

  // insert cohesives
  insertCohesiveElements(facet_nbs_crack_nbs, facet_type, check_only);

  AKANTU_DEBUG_OUT();
  return understressed_candidates.size();
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
  auto & norm_stresses = normal_stresses(type_facet);
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
           make_view(norm_stresses, nb_quad_facet),
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
      auto dist = MeshUtils::distanceBetween2Barycenters(mesh_facets, mesh,
                                                         facet, coh_el);

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
std::tuple<UInt, Real, UInt>
MaterialCohesiveLinearSequential<spatial_dimension>::findCriticalFacetNonLocal(
    const ElementType & type_facet) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model->getMesh();
  const Mesh & mesh_facets = this->model->getMeshFacets();
  // MeshUtils::fillElementToSubElementsData(mesh_facets);
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, spatial_dimension).begin();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  Array<Element> candidate_facets;
  Array<Vector<Real>> final_normals;
  Array<Vector<Real>> facet_normals;
  Array<UInt> coh_crack_nbs;
  Array<Real> highest_stresses;
  Array<UInt> coh_neighbors_nb;
  std::set<Element> crack_contour;
  auto output =
      std::make_tuple(UInt(-1), std::numeric_limits<Real>::min(), UInt(-1));

  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  const auto & facets_check = inserter.getCheckFacets(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  auto & norm_stresses = normal_stresses(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & f_stress = this->model->getStressOnFacets(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & eff_stresses = effective_stresses(type_facet);

  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);

  UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                              .getNbIntegrationPoints(type_cohesive);

  AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                      "The cohesive element and the corresponding facet do "
                      "not have the same numbers of integration points");

  // skip if no facets of this type are present
  UInt nb_facet = f_filter.size();
  if (nb_facet == 0) {
    return output;
  }

  Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
  UInt sp2 = spatial_dimension * spatial_dimension;

  const auto & tangents = this->model->getTangents(type_facet);
  const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                             .getNormalsOnIntegrationPoints(type_facet);
  auto normal_it = normals.begin(spatial_dimension);
  auto tangent_it = tangents.begin(tangents.getNbComponent());
  auto facet_stress_it =
      f_stress.begin(spatial_dimension, spatial_dimension * 2);

  // loop over each facet belonging to this material
  for (auto && data :
       zip(f_filter, sigma_limits, eff_stresses,
           make_view(norm_stresses, nb_quad_facet),
           make_view(norm_tractions, spatial_dimension, nb_quad_facet))) {
    auto facet_nb = std::get<0>(data);
    Element facet{type_facet, facet_nb, _not_ghost};
    auto & sigma_limit = std::get<1>(data);
    auto & eff_stress = std::get<2>(data);
    auto & stress_check = std::get<3>(data);
    auto & normal_traction = std::get<4>(data);

    // skip facets where check shouldn't be inserted (or already inserted)
    if (not facets_check(facet_nb)) {
      continue;
    }

    // compute the effective norm on each quadrature point of the facet
    for (UInt q : arange(nb_quad_facet)) {
      UInt current_quad = facet_nb * nb_quad_facet + q;
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
      stress_check(q) = this->computeEffectiveNorm(stress_tmp, normal, tangent,
                                                   normal_traction_vec);
    }

    // verify if the effective stress overcomes the threshold
    Real final_stress = stress_check.mean();

    // normalize by the limit stress and skip non-stressed facets
    eff_stress = final_stress / sigma_limit;
    if (eff_stress < 1) {
      continue;
    }

    // check to how many cohesive elements this facet is connected
    const Vector<Element> & sub_to_facet =
        mesh_facets.getSubelementToElement(facet);
    std::vector<std::set<Element>> coh_neighbors(sub_to_facet.size());
    std::set<UInt> cohesive_nodes;
    for (auto && data : enumerate(sub_to_facet)) {
      const auto & subel = std::get<1>(data);
      auto & i = std::get<0>(data);
      const auto & connected_elements_dim_min_1 =
          mesh_facets.getElementToSubelement(subel);
      for (const auto & connected_element_dim_min_1 :
           connected_elements_dim_min_1) {
        const auto & connected_elements_dim =
            mesh_facets.getElementToSubelement(connected_element_dim_min_1);
        for (const auto & connected_element_dim : connected_elements_dim) {
          if (mesh.getKind(connected_element_dim.type) == _ek_cohesive) {
            coh_neighbors[i].emplace(connected_element_dim);
            Vector<UInt> coh_nodes =
                mesh.getConnectivity(connected_element_dim);
            for (auto & node : coh_nodes) {
              cohesive_nodes.emplace(node);
            }
            crack_contour.insert(subel);
          }
        }
      }
    }

    Array<UInt> single_tip_subs;
    CSR<Element> nodes_to_cohesives;
    MeshUtils::buildNode2Elements(mesh, nodes_to_cohesives, spatial_dimension,
                                  _ek_cohesive);
    Vector<Real> facet_normal = normal_it[facet_nb * nb_quad_facet];
    Vector<UInt> facet_nodes = mesh_facets.getConnectivity(facet);
    auto coh_normals = this->normals(type_cohesive);
    auto coh_normals_it = coh_normals.begin(spatial_dimension);

    Vector<Real> final_normal(spatial_dimension);
    std::set<Element> visited_cohesives;
    Vector<Real> ref_vector(spatial_dimension);
    for (auto node : cohesive_nodes) {
      for (auto & near_cohesive : nodes_to_cohesives.getRow(node)) {
        if (visited_cohesives.empty()) {
          ref_vector = coh_normals_it[near_cohesive.element * nb_quad_cohesive];
        }
        auto ret = visited_cohesives.emplace(near_cohesive);
        if (ret.second) {
          Vector<Real> coh_normal =
              coh_normals_it[near_cohesive.element * nb_quad_cohesive];
          Real test_dot = ref_vector.dot(coh_normal);
          if (test_dot < 0) {
            coh_normal *= -1;
          }

          final_normal += coh_normal;
        }
      }
    }
    final_normal /= final_normal.norm();

    // compute abs(dot) product between 2 normals & discard sharp angle
    Real dot = facet_normal.dot(final_normal);
    dot = std::abs(dot);

    if (dot < 0.8) {
      goto endloop;
    }
    
    // see if a single tip crack is present
    for (UInt i : arange(sub_to_facet.size())) {
      if (coh_neighbors[i].size() == 1) {
        single_tip_subs.push_back(i);
        // experimental feature (all neighbors have to be single tips)
      } else if (coh_neighbors[i].size() > 1) {
        goto endloop;
      }
    }

    // if no coh els are connected - skip facet
    if (single_tip_subs.empty()) {
      goto endloop;
    }

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
      auto dist = MeshUtils::distanceBetween2Barycenters(mesh_facets, mesh,
                                                         facet, coh_el);

      // ad-hoc rule on barycenters spacing
      // it should discard all elements under sharp angle
      if (dist < 0.8 * facet_indiam) {
          goto endloop;
      }
    }

    // add a candidate facet into a pool
    candidate_facets.push_back(facet);
    coh_crack_nbs.push_back(potential_crack_nb);
    highest_stresses.push_back(eff_stress);
    coh_neighbors_nb.push_back(single_tip_subs.size());
    final_normals.push_back(final_normal);
    facet_normals.push_back(facet_normal);

  endloop:;
  }
  if (candidate_facets.empty()) {
    return output;
  }

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
/*template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::
    insertSingleCohesiveElement(const ElementType & type_facet, UInt facet_nb,
                                UInt facet_crack_nb, bool check_only) {

  Mesh & mesh = this->model->getMesh();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  auto & f_insertion = inserter.getInsertionFacets(type_facet);
  auto & norm_stresses = normal_stresses(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & sig_c_eff = this->sigma_c_eff(type_cohesive);
  auto & del_c = this->delta_c_eff(type_cohesive);
  auto & ins_stress = this->insertion_stress(type_cohesive);
  auto & trac_old = this->tractions.previous(type_cohesive);
  auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type_cohesive);
  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);

  Vector<Real> new_sigmas(nb_quad_facet);
  Array<Real> new_normal_traction(nb_quad_facet, spatial_dimension);
  Vector<Real> new_delta_c(nb_quad_facet);

  // get facet's local id
  auto local_id = f_filter.find(facet_nb);
  AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                      "mismatch between global and local facet numbering"
                          << facet_nb);

  // mark the insertion of the cohesive element
  f_insertion(facet_nb) = true;

  if (check_only)
    return;

  auto & sigma_limit = sigma_limits(local_id);
  Vector<Real> stress_check =
      make_view(norm_stresses, nb_quad_facet).begin()[local_id];
  Matrix<Real> normal_traction =
      make_view(norm_tractions, spatial_dimension, nb_quad_facet)
          .begin()[local_id];

  // store the new cohesive material parameters for each quad point
  for (UInt q = 0; q < nb_quad_facet; ++q) {
    Real new_sigma = stress_check(q);
    Vector<Real> normal_traction_vec(normal_traction(q));

    if (spatial_dimension != 3)
      normal_traction_vec *= -1.;

    new_sigmas(q) = new_sigma;
    for (UInt i : arange(spatial_dimension))
      new_normal_traction(q, i) = normal_traction_vec(i);

    Real new_delta;

    // set delta_c in function of G_c or a given delta_c value
    if (Math::are_float_equal(this->delta_c, 0.))
      new_delta = 2 * this->G_c / new_sigma;
    else
      new_delta = sigma_limit / new_sigma * this->delta_c;

    new_delta_c(q) = new_delta;
  }

  // update material data for the new elements
  UInt old_nb_quad_points = sig_c_eff.size();
  sig_c_eff.resize(old_nb_quad_points + nb_quad_facet);
  ins_stress.resize(old_nb_quad_points + nb_quad_facet);
  trac_old.resize(old_nb_quad_points + nb_quad_facet);
  del_c.resize(old_nb_quad_points + nb_quad_facet);
  crack_numbers.push_back(facet_crack_nb);

  for (UInt q = 0; q < nb_quad_facet; ++q) {
    sig_c_eff(old_nb_quad_points + q) = new_sigmas(q);
    del_c(old_nb_quad_points + q) = new_delta_c(q);
    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      ins_stress(old_nb_quad_points + q, dim) = new_normal_traction(q, dim);
      trac_old(old_nb_quad_points + q, dim) = new_normal_traction(q, dim);
    }
  }
  AKANTU_DEBUG_OUT();
  }*/

/* ----------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::
    insertCohesiveElements(Array<std::pair<UInt, UInt>> & facet_nbs_crack_nbs,
                           ElementType facet_type, bool check_only) {

  Mesh & mesh = this->model->getMesh();
  CohesiveElementInserter & inserter = this->model->getElementInserter();
  ElementType type_cohesive = FEEngine::getCohesiveElementType(facet_type);
  auto & f_filter = this->facet_filter(facet_type);
  auto & f_insertion = inserter.getInsertionFacets(facet_type);
  auto & norm_stresses = normal_stresses(facet_type);
  auto & norm_tractions = normal_tractions(facet_type);
  const auto & sigma_limits = this->sigma_c(facet_type);
  auto & sig_c_eff = this->sigma_c_eff(type_cohesive);
  auto & del_c = this->delta_c_eff(type_cohesive);
  auto & ins_stress = this->insertion_stress(type_cohesive);
  auto & trac_old = this->tractions.previous(type_cohesive);
  auto & crack_numbers = mesh.getData<UInt>("crack_numbers", type_cohesive);
  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(facet_type);
  std::vector<Real> new_sigmas;
  std::vector<Vector<Real>> new_normal_traction;
  std::vector<Real> new_delta_c;

  for (auto && data : facet_nbs_crack_nbs) {
    auto facet_nb = std::get<0>(data);
    auto crack_nb = std::get<1>(data);

    // get facet's local id
    auto local_id = f_filter.find(facet_nb);
    AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                        "mismatch between global and local facet numbering"
                            << facet_nb);

    // mark the insertion of the cohesive element
    f_insertion(facet_nb) = true;

    if (check_only)
      return;

    auto & sigma_limit = sigma_limits(local_id);
    Vector<Real> stress_check =
        make_view(norm_stresses, nb_quad_facet).begin()[local_id];
    Matrix<Real> normal_traction =
        make_view(norm_tractions, spatial_dimension, nb_quad_facet)
            .begin()[local_id];

    // store the new cohesive material parameters for each quad point
    for (UInt q = 0; q < nb_quad_facet; ++q) {
      Real new_sigma = stress_check(q);
      Vector<Real> normal_traction_vec(normal_traction(q));

      if (spatial_dimension != 3)
        normal_traction_vec *= -1.;

      new_sigmas.push_back(new_sigma);
      new_normal_traction.push_back(normal_traction_vec);

      Real new_delta;

      // set delta_c in function of G_c or a given delta_c value
      if (Math::are_float_equal(this->delta_c, 0.))
        new_delta = 2 * this->G_c / new_sigma;
      else
        new_delta = sigma_limit / new_sigma * this->delta_c;

      new_delta_c.push_back(new_delta);
    }

    // insert a corresponding crack number
    crack_numbers.push_back(crack_nb);
  }

  // update material data for the new elements
  UInt old_nb_quad_points = sig_c_eff.size();
  UInt new_nb_quad_points = new_sigmas.size();
  sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
  ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
  trac_old.resize(old_nb_quad_points + new_nb_quad_points);
  del_c.resize(old_nb_quad_points + new_nb_quad_points);

  for (UInt q = 0; q < new_nb_quad_points; ++q) {
    sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
    del_c(old_nb_quad_points + q) = new_delta_c[q];
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
