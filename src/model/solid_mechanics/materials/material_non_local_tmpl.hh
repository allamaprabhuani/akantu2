/**
 * @file   material_non_local_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jul 06 2017
 * @date last modification: Fri Mar 26 2021
 *
 * @brief  Implementation of material non-local
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material.hh"
#include "material_non_local.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, class LocalParent>
MaterialNonLocal<dim, LocalParent>::MaterialNonLocal(
    SolidMechanicsModel & model, const ID & id)
    : LocalParent(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <Int dim, class LocalParent>
void MaterialNonLocal<dim, LocalParent>::insertIntegrationPointsInNeighborhoods(
    GhostType ghost_type,
    const ElementTypeMapReal & quadrature_points_coordinates) {

  IntegrationPoint q;
  q.ghost_type = ghost_type;

  auto & neighborhood = this->model.getNonLocalManager().getNeighborhood(
      this->getNeighborhoodName());

  for (auto & type :
       this->element_filter.elementTypes(dim, ghost_type, _ek_regular)) {
    q.type = type;
    const auto & elem_filter = this->element_filter(type, ghost_type);
    auto nb_element = elem_filter.size();

    if (nb_element == 0) {
      continue;
    }

    auto nb_quad = this->getFEEngine().getNbIntegrationPoints(type, ghost_type);

    auto && quads_it =
        make_view(quadrature_points_coordinates(type, ghost_type), dim, nb_quad)
            .begin();
    for (auto & elem : elem_filter) {
      auto && quads = quads_it[elem];
      q.element = elem;
      for (auto && data : enumerate(quads)) {
        auto nq = std::get<0>(data);
        q.num_point = nq;
        q.global_num = q.element * nb_quad + nq;
        neighborhood.insertIntegrationPoint(q, std::get<1>(data));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, class LocalParent>
void MaterialNonLocal<dim, LocalParent>::updateNonLocalInternals(
    ElementTypeMapReal & non_local_flattened, const ID & field_id,
    GhostType ghost_type, ElementKind kind) {

  /// loop over all types in the material
  for (auto & el_type :
       this->element_filter.elementTypes(dim, ghost_type, kind)) {
    auto & internal =
        this->template getInternal<Real>(field_id)(el_type, ghost_type);

    auto & internal_flat = non_local_flattened(el_type, ghost_type);
    auto nb_component = internal_flat.getNbComponent();

    auto internal_it = internal.begin(nb_component);
    auto internal_flat_it = internal_flat.begin(nb_component);

    /// loop all elements for the given type
    const auto & filter = this->element_filter(el_type, ghost_type);
    Int nb_quads =
        this->getFEEngine().getNbIntegrationPoints(el_type, ghost_type);
    for (auto & elem : filter) {
      for (Int q = 0; q < nb_quads; ++q, ++internal_it) {
        auto global_quad = elem * nb_quads + q;
        *internal_it = internal_flat_it[global_quad];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <Int dim, class LocalParent>
void MaterialNonLocal<dim, LocalParent>::registerNeighborhood() {
  ID name = this->getNeighborhoodName();
  this->model.getNonLocalManager().registerNeighborhood(name, name);
}

} // namespace akantu
