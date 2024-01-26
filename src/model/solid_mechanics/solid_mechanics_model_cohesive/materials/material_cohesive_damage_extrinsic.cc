/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_damage_extrinsic.hh"
#include "dof_synchronizer.hh"
#include "fe_engine_template.hh"
#include "integrator_gauss.hh"
#include "shape_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveDamageExtrinsic<dim>::MaterialCohesiveDamageExtrinsic(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialCohesiveDamage(model, id),
      lambda(registerInternal<Real, CohesiveInternalField>("lambda",
                                                           spatial_dimension)),
      err_openings(registerInternal<Real, CohesiveInternalField>(
          "err_openings", spatial_dimension)),
      czm_damage(
          registerInternal<Real, CohesiveInternalField>("czm_damage", 1))
{
    if(not this->model->getIsExtrinsic())
    {
        AKANTU_EXCEPTION(
            "MaterialCohesiveDamageExtrinsic can be used only with extrinsic cohesive elements");
    }
}

/* -------------------------------------------------------------------------- */
template <Int dim> // NOLINTNEXTLINE(readability-function-cognitive-complexity)
void MaterialCohesiveDamageExtrinsic<dim>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();

//  const auto & mesh_facets = model->getMeshFacets();
//  auto & inserter = model->getElementInserter();

//  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
//    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
//    auto & f_insertion = inserter.getInsertionFacets(type_facet);
//    auto & damage = czm_damage(type_cohesive);

//    const auto & facets_check = inserter.getCheckFacets(type_facet);
//    const auto & facet_filter_array = getFacetFilter(type_facet);

//    auto nb_quad_facet =
//        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

//#ifndef AKANTU_NDEBUG
//    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
//                                .getNbIntegrationPoints(type_cohesive);

//    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
//                        "The cohesive element and the corresponding facet do "
//                        "not have the same numbers of integration points");
//#endif

//    // loop over each facet belonging to this material
//    for (auto && [facet] : facet_filter_array) {
//      // skip facets where check shouldn't be realized
//      if (!facets_check(facet)) {
//        continue;
//      }

//      Vector<Real> d_check(nb_quad_facet);
//      // compute the effective norm on each quadrature point of the facet
//      for (Int q = 0; q < nb_quad_facet; ++q) {
//        auto current_quad = facet * nb_quad_facet + q;
//        auto && normal = normal_begin[current_quad];
//        auto && tangent = tangent_begin[current_quad];
//        auto && facet_stress = facet_stress_begin[current_quad];

//        // compute average stress on the current quadrature point
//        auto && d_1 = facet_stress(0);
//        auto && d_2 = facet_stress(1);

//        auto && stress_avg = (stress_1 + stress_2) / 2.;

//        // compute normal and effective stress
//        stress_check(q) = computeEffectiveNorm(stress_avg, normal, tangent,
//                                               normal_traction(q));
//      }

//      // verify if the effective stress overcomes the threshold
//      auto final_stress = stress_check.mean();
//      if (max_quad_stress_insertion) {
//        final_stress =
//            *std::max_element(stress_check.begin(), stress_check.end());
//      }

//      if (final_stress > sigma_limit) {
//        f_insertion(facet) = true;

//        if (check_only) {
//          continue;
//        }
//      }
//    }
//  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamageExtrinsic<dim>::computeLambdaOnQuad(ElementType type,
                                                      GhostType ghost_type) {
  auto & fem_lambda = model->getFEEngine("LambdaFEEngine");
  const auto & lambda = model->getLambda();
  auto & lambda_on_quad = this->lambda(type, ghost_type);

  auto underlying_type = Mesh::getFacetType(type);
  fem_lambda.interpolateOnIntegrationPoints(
      lambda, lambda_on_quad, dim, underlying_type, ghost_type,
      this->getElementFilter(type, ghost_type));
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamageExtrinsic<dim>::computeTraction(ElementType el_type,
                                                  GhostType ghost_type) {

  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, ghost_type, _ek_cohesive)) {
    auto & elem_filter = getElementFilter(type, ghost_type);
    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    for (auto && args : getArguments(el_type, ghost_type)) {
      computeTractionOnQuad(args);
    }
  }
}

/* -------------------------------------------------------------------------- */
template class MaterialCohesiveDamageExtrinsic<1>;
template class MaterialCohesiveDamageExtrinsic<2>;
template class MaterialCohesiveDamageExtrinsic<3>;
const bool material_is_alocated_cohesive_damage [[maybe_unused]] =
    instantiateMaterial<MaterialCohesiveDamageExtrinsic>("cohesive_damage_extrinsic");

} // namespace akantu
