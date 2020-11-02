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

/* -------------------------------------------------------------------------- */
#include "material_cohesive_linear_sequential.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearSequential<spatial_dimension>::
    MaterialCohesiveLinearSequential(SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<spatial_dimension>(model, id),
      normal_stresses("normal_stresses", *this),
      normal_tractions("normal_tractions", *this),
      effective_stresses("effective_stresses", *this),
      crack_number("crack_number", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam(
      "insertion_threshold", insertion_threshold, Real(0.0),
      _pat_parsable | _pat_readable,
      "Portion of elements below the highest stress to be inserted");

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  normal_stresses.initialize(1);
  normal_tractions.initialize(spatial_dimension);
  effective_stresses.initialize(1);
  crack_number.initialize(1);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::checkInsertion(
    bool check_only) {
  AKANTU_DEBUG_IN();

  Mesh & mesh_facets = this->model->getMeshFacets();
  MeshUtils::fillElementToSubElementsData(mesh_facets);

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
    insertSingleCohesiveElement(type_facet, max_stress_facet, crack_number,
                                check_only);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------*/
template <UInt spatial_dimension>
std::tuple<UInt, Real, Real>
MaterialCohesiveLinearSequential<spatial_dimension>::findCriticalFacet(
    const ElementType & type_facet) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = this->model->getMesh();
  Mesh & mesh_facets = this->model->getMeshFacets();
  MeshUtils::fillElementToSubElementsData(mesh_facets);
  const auto & pos = mesh.getNodes();
  const auto pos_it = make_view(pos, spatial_dimension).begin();
  CohesiveElementInserter & inserter = this->model->getElementInserter();
  CSR<Element> nodes_to_elements;
  MeshUtils::buildNode2Elements(mesh, nodes_to_elements, spatial_dimension,
                                _ek_cohesive);

  UInt max_stress_facet(-1);
  Real max_stress = std::numeric_limits<Real>::min();
  Real potential_crack_nb(-1);

  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
  const auto & facets_check = inserter.getCheckFacets(type_facet);
  auto & f_filter = this->facet_filter(type_facet);
  auto & norm_stresses = normal_stresses(type_facet);
  auto & norm_tractions = normal_tractions(type_facet);
  const auto & f_stress = this->model->getStressOnFacets(type_facet);
  const auto & sigma_limits = this->sigma_c(type_facet);
  auto & eff_stresses = effective_stresses(type_facet);
  auto & facet_conn = mesh_facets.getConnectivity(type_facet);

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
    return std::make_tuple(max_stress_facet, max_stress, potential_crack_nb);

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
    auto facet = std::get<0>(data);
    auto & sigma_limit = std::get<1>(data);
    auto & eff_stress = std::get<2>(data);
    auto & stress_check = std::get<3>(data);
    auto & normal_traction = std::get<4>(data);

    // skip facets where check shouldn't be inserted (or already inserted)
    if (!facets_check(facet))
      continue;

    // check to how many cohesive elements this facet is connected
    std::vector<std::set<Element>> coh_neighbors(2);
    Real coh_crack_nb(-1);
    Array<UInt> facet_nodes(2);
    for (UInt i : arange(2)) {
      facet_nodes(i) = facet_conn(facet, i);
      for (auto & elem : nodes_to_elements.getRow(facet_nodes(i))) {
        coh_neighbors[i].emplace(elem);
      }
    }

    // see if a single tip crack is present
    bool single_tip{false};
    Array<UInt> single_tip_node;
    for (UInt i : arange(2)) {
      if (coh_neighbors[i].size() == 1) {
        single_tip = true;
        single_tip_node.push_back(i);
      }
    }

    // if no coh els are connected or no single tip present - skip facet
    if (not single_tip)
      continue;

    // if single tip is present, compute angle with the candidate facet; if
    // two single tips - compute two angles; if angle(s) are sharp - skip
    bool sharp_angle{false};
    for (auto tip_node : single_tip_node) {
      auto coh_el = *coh_neighbors[tip_node].begin();
      auto & crack_nbs = this->crack_number(coh_el.type, coh_el.ghost_type);
      auto crack_nb_it = crack_nbs.begin();
      auto crack_nb = crack_nb_it[coh_el.element * nb_quad_cohesive];
      coh_crack_nb = crack_nb;

      auto & sub_to_element =
          mesh_facets.getSubelementToElement(coh_el.type, coh_el.ghost_type);
      auto sub_to_el_it = sub_to_element.begin(sub_to_element.getNbComponent());
      const Vector<Element> & subelements_to_element =
          sub_to_el_it[coh_el.element];
      auto & coh_facet_conn =
          mesh_facets.getConnectivity(type_facet, coh_el.ghost_type);
      UInt common_node(-1);
      UInt candidate_facet_node(-1);
      UInt coh_facet_node(-1);
      bool nodes_found{false};

      // iterate through 2 facets of cohesive element
      for (auto & coh_facet : subelements_to_element) {
        Array<UInt> coh_facet_nodes(2);
        for (UInt node : arange(2)) {
          coh_facet_nodes(node) = coh_facet_conn(coh_facet.element, node);
        }
        // find which node is common and which are unique
        for (UInt i : arange(2)) {
          auto id = coh_facet_nodes.find(facet_nodes(i));
          if (id == UInt(-1))
            continue;
          nodes_found = true;
          common_node = facet_nodes(i);
          coh_facet_node = coh_facet_nodes(-1 * id + 1);
          candidate_facet_node = facet_nodes(-1 * i + 1);
        }
        if (nodes_found)
          break;
      }
      AKANTU_DEBUG_ASSERT(nodes_found,
                          "Common and unique nodes were not found between "
                          "candidate facet and a cohesive one");

      // determine the direction of each facet (works only in 2D)
      Vector<Real> coh_facet_dir(spatial_dimension);
      Vector<Real> candidate_facet_dir(spatial_dimension);
      coh_facet_dir = pos_it[coh_facet_node];
      coh_facet_dir -= Vector<Real>(pos_it[common_node]);
      candidate_facet_dir = pos_it[candidate_facet_node];
      candidate_facet_dir -= Vector<Real>(pos_it[common_node]);

      // compute dot product between 2 vectors and discard sharp angles
      Real dot = coh_facet_dir.dot(candidate_facet_dir);
      if (dot >= 0.)
        sharp_angle = true;
    }

    if (sharp_angle)
      continue;

    // compute the effective norm on each quadrature point of the facet
    for (UInt q : arange(nb_quad_facet)) {
      UInt current_quad = facet * nb_quad_facet + q;
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
    if (this->max_quad_stress_insertion)
      final_stress = *std::max_element(stress_check.storage(),
                                       stress_check.storage() + nb_quad_facet);

    // normalize by the limit stress
    eff_stress = final_stress / sigma_limit;

    // update the max stress and corresponding facet
    if (eff_stress > max_stress) {
      max_stress_facet = facet;
      max_stress = eff_stress;
      potential_crack_nb = coh_crack_nb;
    }
  }

  // if (local_max_stress == max_stress)
  return std::make_tuple(max_stress_facet, max_stress, potential_crack_nb);
  // else
  //   return std::make_tuple(-1, -1, -1);
}

/* -------------------------------------------------------------------------*/
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::
    insertSingleCohesiveElement(const ElementType & type_facet, UInt facet_nb,
                                Real crack_nb, bool check_only) {

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
  auto & crack_nbs_current_facet = this->crack_number(type_cohesive);
  UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                           .getNbIntegrationPoints(type_facet);

  Vector<Real> new_sigmas(nb_quad_facet);
  Array<Real> new_normal_traction(nb_quad_facet, spatial_dimension);
  Vector<Real> new_delta_c(nb_quad_facet);
  Vector<Real> new_crack_nb(nb_quad_facet);

  // get facet's local id
  auto local_id = f_filter.find(facet_nb);
  AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                      "mismatch between global and local facet numbering");

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
    new_crack_nb(q) = crack_nb;
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
  crack_nbs_current_facet.resize(old_nb_quad_points + nb_quad_facet);

  for (UInt q = 0; q < nb_quad_facet; ++q) {
    sig_c_eff(old_nb_quad_points + q) = new_sigmas(q);
    del_c(old_nb_quad_points + q) = new_delta_c(q);
    crack_nbs_current_facet(old_nb_quad_points + q) = new_crack_nb(q);
    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      ins_stress(old_nb_quad_points + q, dim) = new_normal_traction(q, dim);
      trac_old(old_nb_quad_points + q, dim) = new_normal_traction(q, dim);
    }
  }
  AKANTU_DEBUG_OUT();
}

INSTANTIATE_MATERIAL(cohesive_linear_sequential,
                     MaterialCohesiveLinearSequential);

} // namespace akantu
