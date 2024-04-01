/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_safe_enum.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
/* -------------------------------------------------------------------------- */

// clang-format off
#ifndef AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_
#define AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Element Types                                                              */
/* -------------------------------------------------------------------------- */
// clang-format off
#define AKANTU_ek_regular_ELEMENT_TYPE  \
  (_point_1)                            \
  (_segment_2)                          \
  (_segment_3)                          \
  (_triangle_3)                         \
  (_triangle_6)                         \
  (_quadrangle_4)                       \
  (_quadrangle_8)                       \
  (_tetrahedron_4)                      \
  (_tetrahedron_10)                     \
  (_pentahedron_6)                      \
  (_pentahedron_15)                     \
  (_hexahedron_8)                       \
  (_hexahedron_20)

#if defined(AKANTU_COHESIVE_ELEMENT)
#define AKANTU_ek_cohesive_ELEMENT_TYPE         \
  (_cohesive_1d_2)                              \
  (_cohesive_2d_4)                              \
  (_cohesive_2d_6)                              \
  (_cohesive_3d_12)                             \
  (_cohesive_3d_16)                             \
  (_cohesive_3d_6)                              \
  (_cohesive_3d_8)
#else
#define AKANTU_ek_cohesive_ELEMENT_TYPE
#endif

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#define AKANTU_ek_structural_ELEMENT_TYPE       \
  (_bernoulli_beam_2)                           \
  (_bernoulli_beam_3)                           \
  (_discrete_kirchhoff_triangle_18)
#else
#define AKANTU_ek_structural_ELEMENT_TYPE
#endif

#if defined(AKANTU_IGFEM)
#define AKANTU_ek_igfem_ELEMENT_TYPE            \
  (_igfem_segment_3)                            \
  (_igfem_triangle_4)                           \
  (_igfem_triangle_5)
#else
#define AKANTU_ek_igfem_ELEMENT_TYPE
#endif

#define AKANTU_ALL_ELEMENT_TYPE                 \
  AKANTU_ek_regular_ELEMENT_TYPE                \
  AKANTU_ek_cohesive_ELEMENT_TYPE               \
  AKANTU_ek_structural_ELEMENT_TYPE             \
  AKANTU_ek_igfem_ELEMENT_TYPE
// clang-format on

/// @enum ElementType type of elements
enum ElementType {
  _not_defined,
  BOOST_PP_SEQ_ENUM(AKANTU_ALL_ELEMENT_TYPE),
  _max_element_type
};

struct ElementType_def {
  using type = ElementType;
  static const type _begin_ = _not_defined;
  static const type _end_ = _max_element_type;
};

using element_type_enum_t = safe_enum<ElementType_def>;
constexpr inline element_type_enum_t element_types{_max_element_type};

/* -------------------------------------------------------------------------- */
/* Element Kinds                                                              */
/* -------------------------------------------------------------------------- */
#define AKANTU_ELEMENT_KIND (_ek_regular)(_ek_cohesive)(_ek_structural)

enum ElementKind { BOOST_PP_SEQ_ENUM(AKANTU_ELEMENT_KIND), _ek_not_defined };

/* -------------------------------------------------------------------------- */
struct ElementKind_def {
  using type = ElementKind;
  static const type _begin_ = BOOST_PP_SEQ_HEAD(AKANTU_ELEMENT_KIND);
  static const type _end_ = _ek_not_defined;
};

using element_kind_t = safe_enum<ElementKind_def>;

/* -------------------------------------------------------------------------- */
/// @enum GeometricalType type of element potentially contained in a Mesh
/// @enum GeometricalType type of element potentially contained in a Mesh
enum GeometricalType {
  _gt_point,
  _gt_segment_2,
  _gt_segment_3,
  _gt_triangle_3,
  _gt_triangle_6,
  _gt_quadrangle_4,
  _gt_quadrangle_8,
  _gt_tetrahedron_4,
  _gt_tetrahedron_10,
  _gt_hexahedron_8,
  _gt_hexahedron_20,
  _gt_pentahedron_6,
  _gt_pentahedron_15,
#if defined(AKANTU_COHESIVE_ELEMENT)
  _gt_cohesive_1d_2,
  _gt_cohesive_2d_4,
  _gt_cohesive_2d_6,
  _gt_cohesive_3d_12,
  _gt_cohesive_3d_16,
  _gt_cohesive_3d_6,
  _gt_cohesive_3d_8,
#endif
#if defined(AKANTU_IGFEM)
  _gt_igfem_segment_3,
  _gt_igfem_triangle_4,
  _gt_igfem_triangle_5,
#endif
  _gt_not_defined
};

/* -------------------------------------------------------------------------- */
/* Interpolation Types                                                        */
/* -------------------------------------------------------------------------- */
// clang-format off
#define AKANTU_itp_regular_INTERPOLATION_TYPES   \
  (_itp_lagrange_point_1)                       \
  (_itp_lagrange_segment_2)                     \
  (_itp_lagrange_segment_3)                     \
  (_itp_lagrange_triangle_3)                    \
  (_itp_lagrange_triangle_6)                    \
  (_itp_lagrange_quadrangle_4)                  \
  (_itp_serendip_quadrangle_8)                  \
  (_itp_lagrange_tetrahedron_4)                 \
  (_itp_lagrange_tetrahedron_10)                \
  (_itp_lagrange_hexahedron_8)                  \
  (_itp_serendip_hexahedron_20)                 \
  (_itp_lagrange_pentahedron_6)                 \
  (_itp_lagrange_pentahedron_15)

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#define AKANTU_itp_structural_INTERPOLATION_TYPES  \
  (_itp_hermite_2)                              \
  (_itp_bernoulli_beam_2)                       \
  (_itp_bernoulli_beam_3)                       \
  (_itp_discrete_kirchhoff_triangle_18)
#else
#define AKANTU_itp_structural_INTERPOLATION_TYPES
#endif

#if defined(AKANTU_IGFEM)
#define AKANTU_itp_igfem_INTERPOLATION_TYPES    \
  (_itp_igfem_segment_3)                        \
  (_itp_igfem_triangle_4)                       \
  (_itp_igfem_triangle_5)
#else
#define AKANTU_itp_igfem_INTERPOLATION_TYPES
#endif

#define AKANTU_INTERPOLATION_TYPES               \
  AKANTU_itp_regular_INTERPOLATION_TYPES         \
  AKANTU_itp_structural_INTERPOLATION_TYPES      \
  AKANTU_itp_igfem_INTERPOLATION_TYPES
// clang-format on

/// @enum InterpolationType type of elements
enum InterpolationType {
  BOOST_PP_SEQ_ENUM(AKANTU_INTERPOLATION_TYPES),
  _itp_not_defined
};

/* -------------------------------------------------------------------------- */
/* Some sub types less probable to change                                     */
/* -------------------------------------------------------------------------- */
/// @enum GeometricalShapeType types of shapes to define the contains
/// function in the element classes
enum GeometricalShapeType {
  _gst_point,
  _gst_triangle,
  _gst_square,
  _gst_prism,
  _gst_not_defined
};

/* -------------------------------------------------------------------------- */
/// @enum GaussIntegrationType classes of types using common
/// description of the gauss point position and weights
enum GaussIntegrationType {
  _git_point,
  _git_segment,
  _git_triangle,
  _git_tetrahedron,
  _git_pentahedron,
  _git_not_defined
};

/* -------------------------------------------------------------------------- */
/// @enum InterpolationKind the family of interpolation types
enum InterpolationKind {
  _itk_lagrangian,
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _itk_structural,
#endif
#if defined(AKANTU_IGFEM)
  _itk_igfem,
#endif
  _itk_not_defined
};

} // namespace akantu

#endif /* AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_ */
// clang-format on

#include "aka_element_classes_info_inline_impl.hh"
