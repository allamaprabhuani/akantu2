/**
 * @file   material_cohesive_bilinear.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Feb 22 2012
 * @date last modification: Sat Dec 19 2020
 *
 * @brief  Bilinear cohesive constitutive law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_bilinear.hh"
//#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
MaterialCohesiveBilinear<spatial_dimension>::MaterialCohesiveBilinear(
    SolidMechanicsModel & model, const ID & id)
    : MaterialCohesiveLinear<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("delta_0", delta_0, Real(0.),
                      _pat_parsable | _pat_readable,
                      "Elastic limit displacement");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  this->sigma_c_eff.setRandomDistribution(this->sigma_c.getRandomParameter());
  MaterialCohesiveLinear<spatial_dimension>::initMaterial();

  this->delta_max.setDefaultValue(delta_0);
  this->insertion_stress.setDefaultValue(0);

  this->delta_max.reset();
  this->insertion_stress.reset();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  AKANTU_DEBUG_IN();

  MaterialCohesiveLinear<spatial_dimension>::onElementsAdded(element_list,
                                                             event);

  bool scale_traction = false;

  // don't scale sigma_c if volume_s hasn't been specified by the user
  if (!Math::are_float_equal(this->volume_s, 0.)) {
    scale_traction = true;
  }

  for (auto && el : element_list) {
    // filter not ghost cohesive elements
    if ((el.ghost_type != _not_ghost) or
        (Mesh::getKind(el.type) != _ek_cohesive)) {
      continue;
    }

    auto index = el.element;
    auto type = el.type;
    auto nb_quad_per_element = this->fem_cohesive.getNbIntegrationPoints(type);

    auto sigma_c_begin =
        make_view(this->sigma_c_eff(type), nb_quad_per_element).begin();
    auto && sigma_c_vec = sigma_c_begin[index];

    auto delta_c_begin =
        make_view(this->delta_c_eff(type), nb_quad_per_element).begin();
    auto && delta_c_vec = delta_c_begin[index];

    if (scale_traction) {
      scaleTraction(el, sigma_c_vec);
    }

    /**
     * Recompute sigma_c as
     * @f$ {\sigma_c}_\textup{new} =
     * \frac{{\sigma_c}_\textup{old} \delta_c} {\delta_c - \delta_0} @f$
     */

    for (Int q = 0; q < nb_quad_per_element; ++q) {
      delta_c_vec(q) = 2 * this->G_c / sigma_c_vec(q);

      if (delta_c_vec(q) - delta_0 < Math::getTolerance()) {
        AKANTU_ERROR("delta_0 = " << delta_0 << " must be lower than delta_c = "
                                  << delta_c_vec(q)
                                  << ", modify your material file");
      }

      sigma_c_vec(q) *= delta_c_vec(q) / (delta_c_vec(q) - delta_0);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
template <typename D1>
void MaterialCohesiveBilinear<spatial_dimension>::scaleTraction(
    const Element & el, Eigen::MatrixBase<D1> & sigma_c_vec) {
  AKANTU_DEBUG_IN();

  auto base_sigma_c = Real(this->sigma_c_eff);

  const auto & mesh_facets = this->model->getMeshFacets();
  const auto & fe_engine = this->model->getFEEngine();

  auto coh_element_to_facet_begin =
      mesh_facets.getSubelementToElement(el.type).begin(2);
  const auto & coh_element_to_facet = coh_element_to_facet_begin[el.element];

  // compute bounding volume
  Real volume = 0.;

  // loop over facets
  for (Int f = 0; f < 2; ++f) {
    const auto & facet = coh_element_to_facet(f);

    const auto & facet_to_element =
        mesh_facets.getElementToSubelement(facet.type, facet.ghost_type);

    // loop over elements connected to each facet
    for (auto && elem : facet_to_element(facet.element)) {
      // skip cohesive elements and dummy elements
      if (elem == ElementNull || Mesh::getKind(elem.type) == _ek_cohesive) {
        continue;
      }

      // unit vector for integration in order to obtain the volume
      auto nb_quadrature_points = fe_engine.getNbIntegrationPoints(elem.type);
      Vector<Real> unit_vector(nb_quadrature_points);
      unit_vector.fill(1);

      volume += fe_engine.integrate(unit_vector, elem);
    }
  }

  // scale sigma_c
  sigma_c_vec = (sigma_c_vec.array() - base_sigma_c) *
                    std::pow(this->volume_s / volume, 1. / this->m_s) +
                base_sigma_c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialCohesiveBilinear<spatial_dimension>::computeTraction(
    const Array<Real> & normal, ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  MaterialCohesiveLinear<spatial_dimension>::computeTraction(normal, el_type,
                                                             ghost_type);

  // adjust damage
  for (auto && data : zip(this->damage(el_type, ghost_type),
                          this->delta_max(el_type, ghost_type),
                          this->delta_c_eff(el_type, ghost_type))) {
    auto & dam = std::get<0>(data);
    auto & delta_max = std::get<1>(data);
    auto & delta_c = std::get<2>(data);

    dam = std::max((delta_max - delta_0) / (delta_c - delta_0), Real(0.));
    dam = std::min(dam, Real(1.));
  }
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(cohesive_bilinear, MaterialCohesiveBilinear);

} // namespace akantu
