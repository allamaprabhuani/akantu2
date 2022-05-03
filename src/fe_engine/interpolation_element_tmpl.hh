/**
 * @file   interpolation_element_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Jun 06 2013
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  interpolation property description
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_INTERPOLATION_ELEMENT_TMPL_HH_
#define AKANTU_INTERPOLATION_ELEMENT_TMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Regular Elements                                                           */
/* -------------------------------------------------------------------------- */
template <> struct InterpolationProperty<_itp_not_defined> {
  static constexpr InterpolationKind kind{_itk_not_defined};
  static constexpr Int nb_nodes_per_element{0};
  static constexpr Int natural_space_dimension{0};
};
template <> struct InterpolationProperty<_itp_lagrange_point_1> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{1};
  static constexpr Int natural_space_dimension{0};
};
template <> struct InterpolationProperty<_itp_lagrange_segment_2> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{2};
  static constexpr Int natural_space_dimension{1};
};
template <> struct InterpolationProperty<_itp_lagrange_segment_3> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{3};
  static constexpr Int natural_space_dimension{1};
};
template <> struct InterpolationProperty<_itp_lagrange_triangle_3> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{3};
  static constexpr Int natural_space_dimension{2};
};
template <> struct InterpolationProperty<_itp_lagrange_triangle_6> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{6};
  static constexpr Int natural_space_dimension{2};
};
template <> struct InterpolationProperty<_itp_lagrange_tetrahedron_4> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{4};
  static constexpr Int natural_space_dimension{3};
};
template <> struct InterpolationProperty<_itp_lagrange_tetrahedron_10> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{10};
  static constexpr Int natural_space_dimension{3};
};
template <> struct InterpolationProperty<_itp_lagrange_quadrangle_4> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{4};
  static constexpr Int natural_space_dimension{2};
};
template <> struct InterpolationProperty<_itp_serendip_quadrangle_8> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{8};
  static constexpr Int natural_space_dimension{2};
};
template <> struct InterpolationProperty<_itp_lagrange_hexahedron_8> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{8};
  static constexpr Int natural_space_dimension{3};
};
template <> struct InterpolationProperty<_itp_serendip_hexahedron_20> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{20};
  static constexpr Int natural_space_dimension{3};
};
template <> struct InterpolationProperty<_itp_lagrange_pentahedron_6> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{6};
  static constexpr Int natural_space_dimension{3};
};
template <> struct InterpolationProperty<_itp_lagrange_pentahedron_15> {
  static constexpr InterpolationKind kind{_itk_lagrangian};
  static constexpr Int nb_nodes_per_element{15};
  static constexpr Int natural_space_dimension{3};
};

} // namespace akantu
#endif /* AKANTU_INTERPOLATION_ELEMENT_TMPL_HH_ */
