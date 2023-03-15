/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <> class GeometricalElement<_gt_igfem_segment_3> {
  static constexpr Int spatial_dimension{1};
  static constexpr UInt nb_nodes_per_element{3};
  static constexpr UInt nb_facet_types{1};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{1}};
  static constexpr std::array<UInt, 2> facet_connectivity_vect{{0, 1}};
};

/* -------------------------------------------------------------------------- */
template <> class GeometricalElement<_gt_igfem_triangle_4> {
  static constexpr Int spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{4};
  static constexpr UInt nb_facet_types{2};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{2, 1}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{2, 3}};
  // clang-format off
  static constexpr std::array<UInt, 7> facet_connectivity_vect{{
      // first type
      0, 2,
      1, 0,
      // second type
      1, 2, 3}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
template <> class GeometricalElement<_gt_igfem_triangle_5> {
  static constexpr Int spatial_dimension{2};
  static constexpr UInt nb_nodes_per_element{5};
  static constexpr UInt nb_facet_types{2};
  static constexpr std::array<UInt, nb_facet_types> nb_facets{{1, 2}};
  static constexpr std::array<UInt, nb_facet_types> nb_nodes_per_facet{{2, 3}};
  // clang-format off
  static constexpr std::array<UInt, 8> facet_connectivity_vect{{
      // first type
      1, 2,
      // second type
      0, 2, 1,
      0, 3, 4}};
  // clang-format on
};

/* -------------------------------------------------------------------------- */
} // namespace akantu
