/**
 * @file   material_cohesive_linear_friction.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Tue Jan 12 2016
 * @date last modification: Wed Feb 21 2018
 *
 * @brief  Linear irreversible cohesive law of mixed mode loading with
 * random stress definition for extrinsic type
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear_friction.hh"
#include "solid_mechanics_model_cohesive.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialCohesiveLinearFriction<spatial_dimension>::
    MaterialCohesiveLinearFriction(SolidMechanicsModel & model, const ID & id)
    : MaterialParent(model, id), mu_insertion("mu_insertion", *this),
      mu_eff("mu_eff", *this), residual_sliding("residual_sliding", *this),
      friction_force("friction_force", *this) {
  AKANTU_DEBUG_IN();

  this->registerParam("mu", mu_insertion, _pat_parsable | _pat_readable,
                      "Value of the friction coefficient at insertion");

  this->registerParam("penalty_for_friction", friction_penalty, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Penalty parameter for the friction behavior");

  this->use_previous_contact_opening = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFriction<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  MaterialParent::initMaterial();

  mu_insertion.initialize(1);
  mu_eff.initialize(1);
  friction_force.initialize(spatial_dimension);
  residual_sliding.initialize(1);
  residual_sliding.initializeHistory();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialCohesiveLinearFriction<spatial_dimension>::checkInsertion(
    bool check_only) {
  AKANTU_DEBUG_IN();

  // get mesh and inserter
  const Mesh & mesh_facets = this->model->getMeshFacets();
  CohesiveElementInserter & inserter = this->model->getElementInserter();

  for (auto && type_facet : mesh_facets.elementTypes(spatial_dimension - 1)) {
    ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);
    const auto & facets_check = inserter.getCheckFacets(type_facet);
    auto & f_insertion = inserter.getInsertionFacets(type_facet);
    auto & f_filter = this->facet_filter(type_facet);
    auto & sig_c_eff = this->sigma_c_eff(type_cohesive);
    auto & mu_eff = this->mu_eff(type_cohesive);
    auto & del_c_eff = this->delta_c_eff(type_cohesive);
    auto & ins_stress = this->insertion_stress(type_cohesive);
    auto & ins_comp = this->insertion_compression(type_cohesive);
    auto & trac_old = this->tractions.previous(type_cohesive);

    UInt nb_quad_facet = this->model->getFEEngine("FacetsFEEngine")
                             .getNbIntegrationPoints(type_facet);

#ifndef AKANTU_NDEBUG
    UInt nb_quad_cohesive = this->model->getFEEngine("CohesiveFEEngine")
                                .getNbIntegrationPoints(type_cohesive);

    AKANTU_DEBUG_ASSERT(nb_quad_cohesive == nb_quad_facet,
                        "The cohesive element and the corresponding facet do "
                        "not have the same numbers of integration points");
#endif

    UInt nb_facet = f_filter.size();
    if (nb_facet == 0)
      continue;

    // get cohesive stress
    const auto & sigma_lim = this->sigma_c(type_facet);
    // get friction at insertion
    const auto & mu_insert = this->mu_insertion(type_facet);
    // get stress on facets
    const auto & f_stress = this->model->getStressOnFacets(type_facet);
    // get normals and tangents to facets
    const auto & tangents = this->model->getTangents(type_facet);
    const auto & normals = this->model->getFEEngine("FacetsFEEngine")
                               .getNormalsOnIntegrationPoints(type_facet);
    // associated iterators
    auto sigma_lim_it = sigma_lim.begin();
    auto mu_insert_it = mu_insert.begin();
    auto facet_stress_begin =
        f_stress.begin(spatial_dimension, spatial_dimension * 2);
    auto normal_begin = normals.begin(spatial_dimension);
    auto tangent_begin = tangents.begin(tangents.getNbComponent());

    // loop utils
    UInt sp2 = spatial_dimension * spatial_dimension;

    Matrix<Real> stress_quad(spatial_dimension, spatial_dimension);
    Matrix<Real> normal_traction_quad(spatial_dimension, nb_quad_facet);
    Vector<Real> normal_compression_quad(nb_quad_facet);
    Vector<Real> effective_stress_quad(nb_quad_facet);
    Vector<Real> mu_insert_quad(nb_quad_facet);
    Vector<Real> sigma_lim_quad(nb_quad_facet);

    std::vector<Real> new_sigmas;
    std::vector<Real> new_mus;
    std::vector<Vector<Real>> new_normal_tractions;
    std::vector<Vector<Real>> new_normal_compressions;
    std::vector<Real> new_delta_cs;

    // loop over each facet belonging to this material
    for (UInt f = 0; f < nb_facet; ++f, ++sigma_lim_it, ++mu_insert_it) {
      // facet id
      UInt facet = f_filter(f);

      // skip facets where check shouldn't be realized
      if (!facets_check(facet))
        continue;

      // initialization for average activation
      Vector<Real> avg_normal(spatial_dimension);
      Vector<Real> avg_tangent(spatial_dimension);
      Vector<Real> avg_normal_traction(spatial_dimension);
      Real avg_normal_compression = 0.;
      Real mu_insert_avg = 0.;
      Real sigma_lim_avg = 0.;

      // compute the effective norm on each quadrature point of the facet
      Real normal_compression_q;
      for (UInt q = 0; q < nb_quad_facet; ++q) {
        UInt current_quad = facet * nb_quad_facet + q;
        const Vector<Real> & normal = normal_begin[current_quad];
        const Vector<Real> & tangent = tangent_begin[current_quad];
        const Matrix<Real> & facet_stress_it = facet_stress_begin[current_quad];

        // compute average stress on the current quadrature point
        Matrix<Real> stress_1(facet_stress_it.storage(), spatial_dimension,
                              spatial_dimension);
        Matrix<Real> stress_2(facet_stress_it.storage() + sp2,
                              spatial_dimension, spatial_dimension);
        stress_quad.copy(stress_1);
        stress_quad += stress_2;
        stress_quad /= 2.;

        // compute normal and effective stress taking into account the sole
        // tensile contribution of the facet tractions
        Vector<Real> normal_traction_vec(normal_traction_quad(q));
        effective_stress_quad(q) = this->computeEffectiveNorm(
            stress_quad, normal, tangent, normal_traction_vec);

        // normal traction on quad
        normal_compression_q = normal_traction_vec.dot(normal);
        // compressive contribution of the normal traction
        normal_compression_quad(q) = std::min(0., normal_compression_q);
        // traction minus compressive part
        normal_traction_quad(q) =
            normal_traction_vec - normal * normal_compression_quad(q);

        // store quad and average infos
        avg_normal_compression += normal_compression_q;
        avg_normal_traction += normal_traction_vec;
        avg_normal += normal;
        avg_tangent += tangent;
        sigma_lim_avg += *sigma_lim_it;
        sigma_lim_quad(q) = *sigma_lim_it;
        mu_insert_avg += *mu_insert_it;
        mu_insert_quad(q) = *mu_insert_it;
      }
      // rescale averages
      sigma_lim_avg /= nb_quad_facet;
      mu_insert_avg /= nb_quad_facet;
      avg_normal /= nb_quad_facet;
      avg_tangent /= nb_quad_facet;
      avg_normal_traction /= nb_quad_facet;
      avg_normal_compression /= nb_quad_facet;

      // compressive contribution of the average normal traction
      avg_normal_compression = std::min(0., avg_normal_compression);
      // average traction minus compressive part
      avg_normal_traction -= avg_normal * avg_normal_compression;

      // verify if the effective stress overcomes the threshold
      Real crit;
      if (this->max_quad_stress_insertion) {
        Vector<Real> crit_quad(nb_quad_facet);
        for (UInt q = 0; q < nb_quad_facet; ++q) {
          crit_quad(q) = effective_stress_quad(q) -
                         (sigma_lim_quad(q) - mu_insert_quad(q) / this->beta *
                                                  normal_compression_quad(q));
        }
        crit = *std::max_element(crit_quad.storage(),
                                 crit_quad.storage() + nb_quad_facet);
        // we must then ensure that the average value we set to
        // new_sigma is above zero
        if (effective_stress_quad.mean() +
                mu_insert_avg / this->beta * avg_normal_compression <
            0) {
          crit = effective_stress_quad.mean() -
                 (sigma_lim_avg -
                  mu_insert_avg / this->beta * avg_normal_compression);
        }
      } else {
        crit = effective_stress_quad.mean() -
               (sigma_lim_avg -
                mu_insert_avg / this->beta * avg_normal_compression);
      }

      if (crit > -Math::getTolerance()) {
        // facet will be inserted
        f_insertion(facet) = true;
        // only if insertion is considered
        if (check_only)
          continue;

        // add tangential friction part to traction
        avg_normal_traction +=
            mu_insert_avg * avg_normal_compression * avg_tangent;
        if (spatial_dimension != 3) {
          avg_normal_traction *= -1.;
          avg_normal_compression *= -1.;
        }

        // store the new cohesive material parameters for each quadrature point
        for (UInt q = 0; q < nb_quad_facet; ++q) {
          // extract tangent
          UInt current_quad = facet * nb_quad_facet + q;
          const Vector<Real> & normal = normal_begin[current_quad];
          const Vector<Real> & tangent = tangent_begin[current_quad];

          // effective stress at insertion
          Real new_sigma =
              effective_stress_quad(q) +
              mu_insert_quad(q) / this->beta * normal_compression_quad(q);

          // friction coefficient
          Real new_mu = mu_insert_quad(q);

          // traction at insertion with friction component
          Vector<Real> new_normal_traction(normal_traction_quad(q));
          new_normal_traction +=
              mu_insert_quad(q) * normal_compression_quad(q) * tangent;
          Vector<Real> new_normal_compression =
              normal_compression_quad(q) * normal;
          if (spatial_dimension != 3) {
            new_normal_traction *= -1.;
            new_normal_compression *= -1.;
          }

          if (new_sigma < 0.05 * sigma_lim_quad(q)) {
            new_mu = mu_insert_avg;
            new_sigma = effective_stress_quad.mean() +
                        mu_insert_avg / this->beta * avg_normal_compression;
            new_normal_traction = avg_normal_traction;
            new_normal_compression = avg_normal_compression * avg_normal;
          }

          // set delta_c in function of G_c or a given delta_c value
          Real new_delta_c;
          if (Math::are_float_equal(this->delta_c, 0.))
            new_delta_c = 2 * this->G_c / new_sigma;
          else
            new_delta_c = (*sigma_lim_it) / new_sigma * this->delta_c;

          // push to insertion list
          new_sigmas.push_back(new_sigma);
          new_mus.push_back(new_mu);
          new_normal_tractions.push_back(new_normal_traction);
          new_normal_compressions.push_back(new_normal_compression);
          new_delta_cs.push_back(new_delta_c);
        }
      }
    }

    // update material data for the new elements
    UInt old_nb_quad_points = sig_c_eff.size();
    UInt new_nb_quad_points = new_sigmas.size();
    sig_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
    mu_eff.resize(old_nb_quad_points + new_nb_quad_points);
    del_c_eff.resize(old_nb_quad_points + new_nb_quad_points);
    ins_stress.resize(old_nb_quad_points + new_nb_quad_points);
    ins_comp.resize(old_nb_quad_points + new_nb_quad_points);
    trac_old.resize(old_nb_quad_points + new_nb_quad_points);

    for (UInt q = 0; q < new_nb_quad_points; ++q) {
      sig_c_eff(old_nb_quad_points + q) = new_sigmas[q];
      mu_eff(old_nb_quad_points + q) = new_mus[q];
      del_c_eff(old_nb_quad_points + q) = new_delta_cs[q];
      for (UInt dim = 0; dim < spatial_dimension; ++dim) {
        ins_stress(old_nb_quad_points + q, dim) = new_normal_tractions[q](dim);
        ins_comp(old_nb_quad_points + q, dim) = new_normal_compressions[q](dim);
        trac_old(old_nb_quad_points + q, dim) = new_normal_tractions[q](dim);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
