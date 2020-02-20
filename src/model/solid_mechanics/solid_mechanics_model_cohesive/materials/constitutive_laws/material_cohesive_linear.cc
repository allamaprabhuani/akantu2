/**
 * @file   material_cohesive_linear.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Thu Jan 14 2021
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear.hh"
#include "dof_synchronizer.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim>
MaterialCohesiveLinear<dim>::MaterialCohesiveLinear(SolidMechanicsModel & model,
                                                    const ID & id)
    : MaterialCohesive(model, id), sigma_c_eff("sigma_c_eff", *this),
      delta_c_eff("delta_c_eff", *this),
      insertion_stress("insertion_stress", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("beta", beta, Real(0.), _pat_parsable | _pat_readable,
                      "Beta parameter");

  this->registerParam("G_c", G_c, Real(0.), _pat_parsable | _pat_readable,
                      "Mode I fracture energy");

  this->registerParam("penalty", penalty, Real(0.),
                      _pat_parsable | _pat_readable, "Penalty coefficient");

  this->registerParam("volume_s", volume_s, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Reference volume for sigma_c scaling");

  this->registerParam("m_s", m_s, Real(1.), _pat_parsable | _pat_readable,
                      "Weibull exponent for sigma_c scaling");

  this->registerParam("kappa", kappa, Real(1.), _pat_parsable | _pat_readable,
                      "Kappa parameter");

  this->registerParam(
      "contact_after_breaking", contact_after_breaking, false,
      _pat_parsable | _pat_readable,
      "Activation of contact when the elements are fully damaged");

  this->registerParam("max_quad_stress_insertion", max_quad_stress_insertion,
                      false, _pat_parsable | _pat_readable,
                      "Insertion of cohesive element when stress is high "
                      "enough just on one quadrature point");

  this->registerParam("recompute", recompute, false,
                      _pat_parsable | _pat_modifiable, "recompute solution");

  this->use_previous_delta_max = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveLinear<dim>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialCohesive::initMaterial();

  sigma_c_eff.initialize(1);
  delta_c_eff.initialize(1);
  insertion_stress.initialize(dim);

  if (not Math::are_float_equal(delta_c, 0.)) {
    delta_c_eff.setDefaultValue(delta_c);
  } else {
    delta_c_eff.setDefaultValue(2 * G_c / sigma_c);
  }

  if (model->getIsExtrinsic()) {
    scaleInsertionTraction();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::updateInternalParameters() {
  /// compute scalars
  beta2_kappa2 = beta * beta / kappa / kappa;
  beta2_kappa = beta * beta / kappa;

  if (Math::are_float_equal(beta, 0)) {
    beta2_inv = 0;
  } else {
    beta2_inv = 1. / beta / beta;
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim> void MaterialCohesiveLinear<dim>::scaleInsertionTraction() {
  AKANTU_DEBUG_IN();

  // do nothing if volume_s hasn't been specified by the user
  if (Math::are_float_equal(volume_s, 0.)) {
    return;
  }

  const auto & mesh_facets = model->getMeshFacets();
  const auto & fe_engine = model->getFEEngine();
  const auto & fe_engine_facet = model->getFEEngine("FacetsFEEngine");

  for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
    const auto & facet_to_element =
        mesh_facets.getElementToSubelement(type_facet);

    auto nb_quad_per_facet = fe_engine_facet.getNbIntegrationPoints(type_facet);

    for(auto && data : enumerate(make_view(sigma_c(type_facet), nb_quad_per_facet))) {
      auto f = std::get<0>(data);
      auto && sigma_c = std::get<1>(data);

      // compute bounding volume
      Real volume = 0;

      for(auto && elem : facet_to_element(f)) {
          if (elem == ElementNull) {
          continue;
        }

        // unit vector for integration in order to obtain the volume
        auto nb_quadrature_points =
            fe_engine.getNbIntegrationPoints(elem.type);
        Vector<Real> unit_vector(nb_quadrature_points);
        unit_vector.fill(1);

        volume += fe_engine.integrate(unit_vector, elem);
      }

      // scale sigma_c
      Vector<Real> base_sigma_c_v(sigma_c.rows());
      sigma_c = (sigma_c.colwise() - base_sigma_c_v) *
                    std::pow(volume_s / volume, 1. / m_s) +
                base_sigma_c_v;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::checkInsertion(bool check_only) {
  AKANTU_DEBUG_IN();

  const auto & mesh_facets = model->getMeshFacets();
  auto & inserter = model->getElementInserter();

  for (auto && type_facet : mesh_facets.elementTypes(dim - 1)) {
    auto type_cohesive = FEEngine::getCohesiveElementType(type_facet);

    const auto & facets_check = inserter.getCheckFacets(type_facet);
    auto & f_insertion = inserter.getInsertionFacets(type_facet);
    auto & f_filter = facet_filter(type_facet);
    auto & sig_c_eff = sigma_c_eff(type_cohesive);
    auto & del_c = delta_c_eff(type_cohesive);
    auto & ins_stress = insertion_stress(type_cohesive);
    auto & trac_old = tractions.previous(type_cohesive);
    const auto & f_stress = model->getStressOnFacets(type_facet);
    const auto & sigma_lim = sigma_c(type_facet);

    auto nb_quad_facet =
        model->getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
    auto nb_quad_cohesive = model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif

    auto nb_facet = f_filter.size();
    //  if (nb_facet == 0) continue;

    auto sigma_lim_it = sigma_lim.begin();

    Matrix<Real, dim, dim> stress_tmp;
    Matrix<Real> normal_traction(dim, nb_quad_facet);
    Vector<Real> stress_check(nb_quad_facet);
    auto sp2 = dim * dim;

    const auto & tangents = model->getTangents(type_facet);
    const auto & normals = model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    auto normal_begin = make_view<dim>(normals).begin();
    auto tangent_begin = make_view<dim, dim == 3 ? 2 : 1>(tangents).begin();
    auto facet_stress_begin = make_view<dim, dim * 2>(f_stress).begin();

    std::vector<Real> new_sigmas;
    std::vector<Vector<Real>> new_normal_traction;
    std::vector<Real> new_delta_c;

    // loop over each facet belonging to this material
    for (Int f = 0; f < nb_facet; ++f, ++sigma_lim_it) {
      auto facet = f_filter(f);
      // skip facets where check shouldn't be realized
      if (!facets_check(facet)) {
        continue;
      }

      // compute the effective norm on each quadrature point of the facet
      for (Int q = 0; q < nb_quad_facet; ++q) {
        auto current_quad = facet * nb_quad_facet + q;
        auto && normal = normal_begin[current_quad];
        auto && tangent = tangent_begin[current_quad];
        auto && facet_stress_it = facet_stress_begin[current_quad];

        // compute average stress on the current quadrature point
        MatrixProxy<const Real, dim, dim> stress_1(facet_stress_it.data());
        MatrixProxy<const Real, dim, dim> stress_2(facet_stress_it.data() +
                                                   sp2);

        auto && stress_avg = (stress_1 + stress_2) / 2.;

        // compute normal and effective stress
        stress_check(q) = computeEffectiveNorm(stress_avg, normal, tangent,
                                               normal_traction(q));
      }

      // verify if the effective stress overcomes the threshold
      auto final_stress = stress_check.mean();
      if (max_quad_stress_insertion) {
        final_stress = *std::max_element(
            stress_check.data(), stress_check.data() + nb_quad_facet);
      }

      if (final_stress > *sigma_lim_it) {
        f_insertion(facet) = true;

        if (check_only) {
          continue;
        }

        // store the new cohesive material parameters for each quadrature
        // point
        for (Int q = 0; q < nb_quad_facet; ++q) {
          auto new_sigma = stress_check(q);
          auto && normal_traction_vec = normal_traction(q);

          if (dim != 3) {
            normal_traction_vec *= -1.;
          }

          new_sigmas.push_back(new_sigma);
          new_normal_traction.push_back(normal_traction_vec);

          Real new_delta;

          // set delta_c in function of G_c or a given delta_c value
          if (Math::are_float_equal(delta_c, 0.)) {
            new_delta = 2 * G_c / new_sigma;
          } else {
            new_delta = (*sigma_lim_it) / new_sigma * delta_c;
          }

          new_delta_c.push_back(new_delta);
        }
      }
    }

    // update material data for the new elements
    auto old_nb_quad_points = sig_c_eff.size();
    auto new_nb_quad_points = new_sigmas.size();
    sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
    ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
    trac_old.resize(old_nb_quad_points + new_nb_quad_points);
    del_c.resize(old_nb_quad_points + new_nb_quad_points);

    for (UInt q = 0; q < new_nb_quad_points; ++q) {
      sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
      del_c(old_nb_quad_points + q) = new_delta_c[q];
      for (Int d = 0; d < dim; ++d) {
        ins_stress(old_nb_quad_points + q, d) = new_normal_traction[q](d);
        trac_old(old_nb_quad_points + q, d) = new_normal_traction[q](d);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::computeTraction(const Array<Real> & normal,
                                                  ElementType el_type,
                                                  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  auto traction_it = tractions(el_type, ghost_type).begin(dim);
  auto opening_it = opening(el_type, ghost_type).begin(dim);
  auto contact_traction_it = contact_tractions(el_type, ghost_type).begin(dim);
  auto contact_opening_it = contact_opening(el_type, ghost_type).begin(dim);

  auto normal_it = normal.begin(dim);
  auto traction_end = tractions(el_type, ghost_type).end(dim);
  auto sigma_c_it = sigma_c_eff(el_type, ghost_type).begin();
  auto delta_max_it = delta_max(el_type, ghost_type).begin();
  auto delta_c_it = delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = damage(el_type, ghost_type).begin();
  auto insertion_stress_it = insertion_stress(el_type, ghost_type).begin(dim);

  Vector<Real> normal_opening(dim);
  Vector<Real> tangential_opening(dim);

  /// loop on each quadrature point
  for (; traction_it != traction_end;
       ++traction_it, ++opening_it, ++normal_it, ++sigma_c_it, ++delta_max_it,
       ++delta_c_it, ++damage_it, ++contact_traction_it, ++insertion_stress_it,
       ++contact_opening_it) {
    Real normal_opening_norm{0};
    Real tangential_opening_norm{0};
    bool penetration{false};
    this->computeTractionOnQuad(
        *traction_it, *opening_it, *normal_it, *delta_max_it, *delta_c_it,
        *insertion_stress_it, *sigma_c_it, normal_opening, tangential_opening,
        normal_opening_norm, tangential_opening_norm, *damage_it, penetration,
        *contact_traction_it, *contact_opening_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim>
void MaterialCohesiveLinear<dim>::computeTangentTraction(
    ElementType el_type, Array<Real> & tangent_matrix,
    const Array<Real> & normal, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /// define iterators
  auto tangent_it = tangent_matrix.begin(dim, dim);

  auto tangent_end = tangent_matrix.end(dim, dim);

  auto normal_it = normal.begin(dim);

  auto opening_it = opening(el_type, ghost_type).begin(dim);

  /// NB: delta_max_it points on delta_max_previous, i.e. the
  /// delta_max related to the solution of the previous incremental
  /// step
  auto delta_max_it = delta_max.previous(el_type, ghost_type).begin();
  auto sigma_c_it = sigma_c_eff(el_type, ghost_type).begin();
  auto delta_c_it = delta_c_eff(el_type, ghost_type).begin();
  auto damage_it = damage(el_type, ghost_type).begin();
  auto contact_opening_it = contact_opening(el_type, ghost_type).begin(dim);

  Vector<Real> normal_opening(dim);
  Vector<Real> tangential_opening(dim);

  for (; tangent_it != tangent_end; ++tangent_it, ++normal_it, ++opening_it,
                                    ++delta_max_it, ++sigma_c_it, ++delta_c_it,
                                    ++damage_it, ++contact_opening_it) {
    Real normal_opening_norm{0};
    Real tangential_opening_norm{0};
    bool penetration{false};
    this->computeTangentTractionOnQuad(
        *tangent_it, *delta_max_it, *delta_c_it, *sigma_c_it, *opening_it,
        *normal_it, normal_opening, tangential_opening, normal_opening_norm,
        tangential_opening_norm, *damage_it, penetration, *contact_opening_it);
  }

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_linear, MaterialCohesiveLinear);

} // namespace akantu
