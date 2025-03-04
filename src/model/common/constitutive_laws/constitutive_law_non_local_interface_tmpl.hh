/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "constitutive_law_non_local_interface.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int dim, class ConstitutiveLawNonLocalInterface,
          class ConstitutiveLawParent>
ConstitutiveLawNonLocal<dim, ConstitutiveLawNonLocalInterface,
                        ConstitutiveLawParent>::
    ConstitutiveLawNonLocal(
        typename ConstitutiveLawParent::ConstitutiveLawsHandler & handler,
        const ID & id)
    : ConstitutiveLawParent(handler, id) {}

/* -------------------------------------------------------------------------- */
template <Int dim, class ConstitutiveLawNonLocalInterface,
          class ConstitutiveLawParent>
void ConstitutiveLawNonLocal<dim, ConstitutiveLawNonLocalInterface,
                             ConstitutiveLawParent>::
    insertIntegrationPointsInNeighborhoods(
        GhostType ghost_type,
        const ElementTypeMapReal & quadrature_points_coordinates) {

  IntegrationPoint q;
  q.ghost_type = ghost_type;

  auto & neighborhood = this->getModel().getNonLocalManager().getNeighborhood(
      this->getNeighborhoodName());

  for (auto && type :
       this->getElementFilter().elementTypes(dim, ghost_type, _ek_regular)) {
    q.type = type;
    const auto & elem_filter = this->getElementFilter(type, ghost_type);
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
template <Int dim, class ConstitutiveLawNonLocalInterface,
          class ConstitutiveLawParent>
void ConstitutiveLawNonLocal<dim, ConstitutiveLawNonLocalInterface,
                             ConstitutiveLawParent>::
    updateNonLocalInternals(ElementTypeMapReal & non_local_flattened,
                            const ID & field_id, GhostType ghost_type,
                            ElementKind kind) {

  /// loop over all types in the material
  for (auto && el_type :
       this->getElementFilter().elementTypes(dim, ghost_type, kind)) {
    auto & internal =
        this->template getInternal<Real>(field_id)(el_type, ghost_type);

    auto & internal_flat = non_local_flattened(el_type, ghost_type);
    auto nb_component = internal_flat.getNbComponent();

    auto internal_it = internal.begin(nb_component);
    auto internal_flat_it = internal_flat.begin(nb_component);

    /// loop all elements for the given type
    const auto & filter = this->getElementFilter(el_type, ghost_type);
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
template <Int dim, class ConstitutiveLawNonLocalInterface,
          class ConstitutiveLawParent>
void ConstitutiveLawNonLocal<dim, ConstitutiveLawNonLocalInterface,
                             ConstitutiveLawParent>::registerNeighborhood() {
  ID name = this->getNeighborhoodName();
  this->handler.getNonLocalManager().registerNeighborhood(name, name);
}

} // namespace akantu
