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
      lambda(registerInternal<Real, FacetInternalField>("lambda",
                                                           spatial_dimension)),
      czm_damage(
          registerInternal<Real, FacetInternalField>("czm_damage", 1)),
      insertion_stress(registerInternal<Real, CohesiveInternalField>(
          "insertion_stress", dim)),
      is_new_crack(registerInternal<Real, CohesiveInternalField>(
          "is_new_crack", 1))
{
    if(not this->model->getIsExtrinsic())
    {
        AKANTU_EXCEPTION(
            "MaterialCohesiveDamageExtrinsic can be used only with extrinsic cohesive elements");
    }
//    czm_damage.setDefaultValue(0.1);
    is_new_crack.setDefaultValue(0.);
}

/* -------------------------------------------------------------------------- */
template <Int dim> // NOLINTNEXTLINE(readability-function-cognitive-complexity)
void MaterialCohesiveDamageExtrinsic<dim>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();

  const auto & mesh_facets = model->getMeshFacets();
  auto & inserter = model->getElementInserter();

  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    auto & f_insertion = inserter.getInsertionFacets(type_facet);
    auto & lambda = this->lambda(type_facet);
    auto & damage = this->czm_damage(type_facet);
    auto & ins_stress = insertion_stress(type_cohesive);
    auto & is_new_crack = this->is_new_crack(type_cohesive);
    Int old_nb_quad_points = is_new_crack.size();
    Int new_nb_quad_points = 0;

    const auto & facets_check = inserter.getCheckFacets(type_facet);
    const auto & facet_filter_array = getFacetFilter(type_facet);

    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif

    auto damage_begin = damage.begin();
    auto lambda_begin = make_view<dim>(lambda).begin();

    std::vector<Vector<Real>> new_normal_traction;

    // loop over each facet belonging to this material
    for (auto && facet : facet_filter_array) {
      // skip facets where check shouldn't be realized
      if (!facets_check(facet)) {
        continue;
      }

      // compute the effective norm on each quadrature point of the facet
      Real d(0.);
      for (Int q = 0; q < nb_quad_facet; ++q) {
        auto current_quad = facet * nb_quad_facet + q;
        d+=damage_begin[current_quad];
      }
      d/=nb_quad_facet;

      // verify if the effective stress overcomes the threshold
      if (not Math::are_float_equal(d, 0.)) {
        f_insertion(facet) = true;
        if (check_only) {
          continue;
        }
        for (Int q = 0; q < nb_quad_facet; ++q) {
           auto current_quad = facet * nb_quad_facet + q;
           new_normal_traction.emplace_back(lambda_begin[current_quad]);
           new_nb_quad_points++;
        }
      }
    }
    // update material data for the new elements
    ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
    is_new_crack.resize(old_nb_quad_points + new_nb_quad_points);
    for (Int q = 0; q < new_nb_quad_points; ++q) {
      is_new_crack(old_nb_quad_points + q) = 1;
      for (Int d = 0; d < dim; ++d) {
        ins_stress(old_nb_quad_points + q, d) = new_normal_traction[q](d);
      }
    }
  }

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

    std::cout << " In MaterialCohesiveDamageExtrinsic<dim>::computeTraction " << std::endl ;

  const auto & mesh_facets = model->getMeshFacets();
  for (const auto & type : getElementFilter().elementTypes(
           spatial_dimension, ghost_type, _ek_cohesive)) {
    auto type_facet = Mesh::getFacetType(type);

    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif

    auto & elem_filter = this->getElementFilter(type, ghost_type);

    auto nb_element = elem_filter.size();
    if (nb_element == 0) {
      continue;
    }

    const auto & element_to_facet =
        mesh_facets.getSubelementToElement(el_type,ghost_type);

    auto insertion_stress_begin = make_view<dim>(insertion_stress(el_type,ghost_type)).begin();
    auto is_new_crack_begin = is_new_crack(el_type,ghost_type).begin();

    auto damage_begin = czm_damage(type_facet,ghost_type).begin();
    auto opening_begin = make_view<dim>(opening(el_type,ghost_type)).begin();
    auto traction_begin = make_view<dim>(tractions(el_type,ghost_type)).begin();

    for (auto && elem : elem_filter) {

      auto && facet = -1;
      auto && facet_1 = element_to_facet(elem,0);
      auto && facet_2 = element_to_facet(elem,1);
      facet = facet_1<facet_2?facet_1.element:facet_2.element;
      for (Int q = 0; q < nb_quad_facet; ++q) {
        auto current_quad = elem * nb_quad_facet + q;
        auto && is_new = is_new_crack_begin[current_quad];
        auto && t = traction_begin[current_quad];
        if(is_new > 0)
        {
          t = insertion_stress_begin[current_quad];
          is_new = 0;
        }
        else
        {
          auto current_facet_quad = facet * nb_quad_facet + q;
          auto && w = opening_begin[current_quad];
          auto && d = damage_begin[current_facet_quad];
          std::cout << "facet = " << facet << std::endl ;
          std::cout << " w = " << w << std::endl ;
          std::cout << " d = " << d << std::endl;
          t = k*stiffness(d)*w;
          std::cout << " t = " << t << std::endl;
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveDamageExtrinsic<dim>::computeLambda(GhostType ghost_type)
{
  std::cout << " In MaterialCohesiveDamageExtrinsic<dim>::computeLambda " << std::endl ;
  const auto & mesh_facets = model->getMeshFacets();
  for (const auto & type_facet : mesh_facets.elementTypes(dim - 1)) {
    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

    #ifndef AKANTU_NDEBUG
    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
    #endif

    const auto & element_to_subelement =
            mesh_facets.getElementToSubelement(type_facet, ghost_type);

    const auto & normals = model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    const auto & f_stress = model->getStressOnFacets(type_facet);

    auto facet_stress_begin = make_view(f_stress, dim, dim, 2).begin();
    auto normal_begin = make_view<dim>(normals).begin();

    auto damage_begin = czm_damage(type_facet,ghost_type).begin();
    auto opening_begin = make_view<dim>(opening(type_cohesive,ghost_type)).begin();

    for (auto data :
         enumerate(make_view(lambda(type_facet), dim,nb_quad_facet))) {
      auto && f = std::get<0>(data);
      auto && lda = std::get<1>(data);
      auto && elem = element_to_subelement(f);

      auto && e = -1;
      for(auto ec : elem)
      {
        if(ec!=ElementNull)
        {
          if(ec.kind() == _ek_cohesive)
          {
              e = ec.element;
              continue;
          }
        }
      }

      if(e < 0)
      {
        for (Int q = 0; q < nb_quad_facet; ++q) {
          auto current_quad = f * nb_quad_facet + q;
          auto && normal = normal_begin[current_quad];
          auto && facet_stress = facet_stress_begin[current_quad];

          auto && stress_1 = facet_stress(0);
          auto && stress_2 = facet_stress(1);

          auto && stress_avg = (stress_1 + stress_2) / 2.;
          lda(q) = stress_avg*normal;
        }
      }
      else
      {
        for (Int q = 0; q < nb_quad_facet; ++q) {
          auto current_quad = f * nb_quad_facet + q;
          auto current_elem_quad = e * nb_quad_facet + q;
          auto && w = opening_begin[current_elem_quad];
          auto && d = damage_begin[current_quad];
          lda(q) = k*augmented_stiffness(d)*w;
          std::cout << "f = " << f << std::endl ;
          std::cout << " w = " << w << std::endl ;
          std::cout << " d = " << d << std::endl ;
          std::cout << " lda = " << lda(q) << std::endl ;
        }
      }
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
