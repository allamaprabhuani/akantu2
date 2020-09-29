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
      effective_stresses("effective_stresses", *this) {
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
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearSequential<spatial_dimension>::checkInsertion(
    bool check_only) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh_facets = this->model->getMeshFacets();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  for (auto && type_facet : mesh_facets.elementTypes(spatial_dimension - 1)) {
    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    const auto & facets_check = inserter.getCheckFacets(type_facet);
    auto & f_insertion = inserter.getInsertionFacets(type_facet);
    auto & f_filter = this->facet_filter(type_facet);
    auto & sig_c_eff = this->sigma_c_eff(type_cohesive);
    auto & del_c = this->delta_c_eff(type_cohesive);
    auto & ins_stress = this->insertion_stress(type_cohesive);
    auto & trac_old = this->tractions.previous(type_cohesive);
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

    Matrix<Real> stress_tmp(spatial_dimension, spatial_dimension);
    UInt sp2 = spatial_dimension * spatial_dimension;

    const auto & tangents = this->model->getTangents(type_facet);
    const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    auto normal_begin = normals.begin(spatial_dimension);
    auto tangent_begin = tangents.begin(tangents.getNbComponent());
    auto facet_stress_begin =
        f_stress.begin(spatial_dimension, spatial_dimension * 2);

    std::vector<Real> new_sigmas;
    std::vector<Vector<Real>> new_normal_traction;
    std::vector<Real> new_delta_c;
    std::vector<std::pair<UInt, Real>> facet_stress;

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

      // skip facets where check shouldn't be realized
      if (!facets_check(facet))
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
      if (this->max_quad_stress_insertion)
        final_stress = *std::max_element(
            stress_check.storage(), stress_check.storage() + nb_quad_facet);

      // normalize by the limit stress
      eff_stress = final_stress / sigma_limit;
      facet_stress.push_back(std::make_pair(facet, eff_stress));
    }

    // sort facets by the value of effective stress
    std::sort(facet_stress.begin(), facet_stress.end(),
              [](const std::pair<UInt, Real> & lhs,
                 const std::pair<UInt, Real> & rhs) {
                return lhs.second > rhs.second;
              });

    // continue to the next el type if no elements to insert
    if (facet_stress[0].second <= 1.)
      continue;

    // second loop to activate certain portion of elements
    const Mesh & mesh = this->model->getMesh();
    CSR<UInt> nodes_to_elements;
    MeshUtils::buildNode2ElementsElementTypeMap(mesh, nodes_to_elements,
                                                _cohesive_2d_4, _not_ghost);

    // for (auto && data :
    //      zip(f_filter, sigma_limits, eff_stresses,
    //          make_view(norm_stresses, nb_quad_facet),
    //          make_view(norm_tractions, spatial_dimension, nb_quad_facet))) {
    UInt nb_coh_inserted{0};
    for (auto && pair : facet_stress) {
      // allow insertion of only 1 cohesive element
      if (nb_coh_inserted == 1)
        break;

      auto & facet = std::get<0>(pair);
      // get facet's local id
      auto local_id = f_filter.find(facet);
      AKANTU_DEBUG_ASSERT(local_id != UInt(-1),
                          "mismatch between global and local facet numbering");

      // skip facets where check shouldn't be realized
      if (!facets_check(facet))
        continue;

      auto & eff_stress = std::get<1>(pair);
      // break when descending eff_stresses drop below 1
      if (eff_stress < 1.)
        break;

      auto & sigma_limit = sigma_limits(local_id);
      Vector<Real> stress_check =
          make_view(norm_stresses, nb_quad_facet).begin()[local_id];
      Matrix<Real> normal_traction =
          make_view(norm_tractions, spatial_dimension, nb_quad_facet)
              .begin()[local_id];

      // Real reduced_threshold =
      //     max_eff_stress - (max_eff_stress - 1) * this->insertion_threshold;
      // reduced_threshold = std::max(reduced_threshold, 1.);

      // if (eff_stress >= reduced_threshold) {
      // check if this facet is connected to a single cohesive else
      std::set<UInt> coh_neighbors;
      for (UInt i : arange(2)) {
        const UInt facet_node = facet_conn(facet, i);
        for (auto & elem : nodes_to_elements.getRow(facet_node)) {
          // if (mesh.getKind(elem.type) == _ek_cohesive)
          coh_neighbors.emplace(elem);
        }
      }
      // if no or more than one coh neighbors - continue
      if (coh_neighbors.size() != 1)
        continue;

      f_insertion(facet) = true;
      nb_coh_inserted++;

      if (check_only)
        continue;

      // store the new cohesive material parameters for each quadrature
      // point
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
      // }
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
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_linear_sequential,
                     MaterialCohesiveLinearSequential);

} // namespace akantu
