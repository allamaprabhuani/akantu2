/**
 * @file   element_class.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Declaration of the ElementClass main class and the
 * Integration and Interpolation elements
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
#include "aka_common.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_HH_
#define AKANTU_ELEMENT_CLASS_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/// default element class structure
template <ElementType element_type> struct ElementClassProperty {
  static constexpr GeometricalType geometrical_type{_gt_not_defined};
  static constexpr InterpolationType interpolation_type{_itp_not_defined};
  static constexpr ElementKind element_kind{_ek_regular};
  static constexpr Int spatial_dimension{0};
  static constexpr GaussIntegrationType gauss_itegration_type{_git_not_defined};
  static constexpr Int polynomial_degree{0};
};

#if !defined(DOXYGEN)
/// Macro to generate the element class structures for different element types
#define AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(elem_type, geom_type,             \
                                             interp_type, elem_kind, sp,       \
                                             gauss_int_type, min_int_order)    \
  template <> struct ElementClassProperty<elem_type> {                         \
    static constexpr GeometricalType geometrical_type{geom_type};              \
    static constexpr InterpolationType interpolation_type{interp_type};        \
    static constexpr ElementKind element_kind{elem_kind};                      \
    static constexpr Int spatial_dimension{sp};                                \
    static constexpr GaussIntegrationType gauss_integration_type{              \
        gauss_int_type};                                                       \
    static constexpr Int polynomial_degree{min_int_order};                     \
  }
#else
#define AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(elem_type, geom_type,             \
                                             interp_type, elem_kind, sp,       \
                                             gauss_int_type, min_int_order)
#endif

/* -------------------------------------------------------------------------- */
/* Geometry                                                                   */
/* -------------------------------------------------------------------------- */
/// Default GeometricalShape structure
template <GeometricalType geometrical_type> struct GeometricalShape {
  static constexpr GeometricalShapeType shape{_gst_point};
};

/// Templated GeometricalShape with function contains
template <GeometricalShapeType shape> struct GeometricalShapeContains {
  /// Check if the point (vector in 2 and 3D) at natural coordinate coord
  template <class D>
  static inline bool contains(const Eigen::MatrixBase<D> &coord);
};

#if !defined(DOXYGEN)
/// Macro to generate the GeometricalShape structures for different geometrical
/// types
#define AKANTU_DEFINE_SHAPE(geom_type, geom_shape)                             \
  template <> struct GeometricalShape<geom_type> {                             \
    static constexpr GeometricalShapeType shape{geom_shape};                   \
  }

AKANTU_DEFINE_SHAPE(_gt_hexahedron_20, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_hexahedron_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_pentahedron_15, _gst_prism);
AKANTU_DEFINE_SHAPE(_gt_pentahedron_6, _gst_prism);
AKANTU_DEFINE_SHAPE(_gt_point, _gst_point);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_4, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_quadrangle_8, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_segment_2, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_segment_3, _gst_square);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_10, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_tetrahedron_4, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_triangle_3, _gst_triangle);
AKANTU_DEFINE_SHAPE(_gt_triangle_6, _gst_triangle);
#endif
/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type>
struct GeometricalElementProperty {};

template <ElementType element_type>
struct ElementClassExtraGeometryProperties {};

/* -------------------------------------------------------------------------- */
/// Templated GeometricalElement with function getInradius
template <GeometricalType geometrical_type,
          GeometricalShapeType shape =
              GeometricalShape<geometrical_type>::shape>
class GeometricalElement {
  using geometrical_property = GeometricalElementProperty<geometrical_type>;

public:
  /// compute the in-radius: \todo should be renamed for characteristic length
  template <class D>
  static inline Real getInradius(const Eigen::MatrixBase<D> & /*X*/) {
    return 0.;
  }

  /// compute the normal to the element
  template <class D1, class D2>
  static inline void getNormal(const Eigen::MatrixBase<D1> & /*X*/,
                               Eigen::MatrixBase<D2> &n) {
    n.zero();
  }

  /// true if the natural coordinates are in the element
  template <class D>
  static inline bool contains(const Eigen::MatrixBase<D> &coord);

public:
  static constexpr auto getSpatialDimension() {
    return geometrical_property::spatial_dimension;
  }
  static constexpr auto getNbNodesPerElement() {
    return geometrical_property::nb_nodes_per_element;
  }
  static inline constexpr auto getNbFacetTypes() {
    return geometrical_property::nb_facet_types;
  };
  static inline constexpr Int getNbFacetsPerElement(Idx t);
  static inline constexpr Int getNbFacetsPerElement();

  static inline constexpr decltype(auto)
  getFacetLocalConnectivityPerElement(Idx t = 0);

  template <Idx t,
            std::size_t size = std::tuple_size<
                decltype(geometrical_property::nb_facets)>::value,
            std::enable_if_t<(t < size)> * = nullptr>
  static inline constexpr decltype(auto) getFacetLocalConnectivityPerElement();

  template <Idx t,
            std::size_t size = std::tuple_size<
                decltype(geometrical_property::nb_facets)>::value,
            std::enable_if_t<not(t < size)> * = nullptr>
  static inline constexpr decltype(auto) getFacetLocalConnectivityPerElement();
};

/* -------------------------------------------------------------------------- */
/* Interpolation                                                              */
/* -------------------------------------------------------------------------- */
/// default InterpolationProperty structure
template <InterpolationType interpolation_type> struct InterpolationProperty {};

#if !defined(DOXYGEN)
/// Macro to generate the InterpolationProperty structures for different
/// interpolation types
#define AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(itp_type, itp_kind,          \
                                                  nb_nodes, ndim)              \
  template <> struct InterpolationProperty<itp_type> {                         \
    static constexpr InterpolationKind kind{itp_kind};                         \
    static constexpr Int nb_nodes_per_element{nb_nodes};                       \
    static constexpr Int natural_space_dimension{ndim};                        \
  }
#else
#define AKANTU_DEFINE_INTERPOLATION_TYPE_PROPERTY(itp_type, itp_kind,          \
                                                  nb_nodes, ndim)
#endif
/* -------------------------------------------------------------------------- */
/// Generic (templated by the enum InterpolationType which specifies the order
/// and the dimension of the interpolation) class handling the elemental
/// interpolation
template <InterpolationType interpolation_type,
          InterpolationKind kind =
              InterpolationProperty<interpolation_type>::kind>
class InterpolationElement {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;

  /// compute the shape values for a given set of points in natural coordinates
  template <class D1, class D2,
            aka::enable_if_t<aka::are_matrices<D1, D2>::value> * = nullptr>
  static inline void computeShapes(const Eigen::MatrixBase<D1> &Xs,
                                   const Eigen::MatrixBase<D2> &N_);

  /// compute the shape values for a given point in natural coordinates
  template <class D1, class D2,
            aka::enable_if_t<aka::are_vectors<D1, D2>::value> * = nullptr>
  static inline void computeShapes(const Eigen::MatrixBase<D1> & /*Xs*/,
                                   Eigen::MatrixBase<D2> & /*N_*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  template <class D>
  static inline void computeDNDS(const Eigen::MatrixBase<D> &Xs,
                                 Tensor3Base<Real> &dnds);
  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  template <class D1, class D2>
  static inline void computeDNDS(const Eigen::MatrixBase<D1> & /*Xs*/,
                                 Eigen::MatrixBase<D2> & /*dNdS*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /**
   * compute @f$ @f$

   **/
  static inline void computeD2NDS2(const Matrix<Real> &natural_coord,
                                   Tensor3<Real> &d2nds2);

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the
   * second variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  template <class vector_type, class matrix_type>
  static inline void computeD2NDS2(const vector_type & /*unused*/,
                                   matrix_type & /*unused*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  /// in the case of spatial_dimension != natural_space_dimension
  template <class D>
  static inline Real computeSpecialJacobian(const Eigen::MatrixBase<D> &) {
    AKANTU_TO_IMPLEMENT();
  }

  /// interpolate a field given (arbitrary) natural coordinates
  template <class Derived1, class Derived2>
  static inline decltype(auto) interpolateOnNaturalCoordinates(
      const Eigen::MatrixBase<Derived1> &natural_coords,
      const Eigen::MatrixBase<Derived2> &nodal_values) {
    using interpolation = InterpolationProperty<interpolation_type>;
    Eigen::Matrix<Real, interpolation::nb_nodes_per_element, 1> shapes;
    computeShapes(natural_coords, shapes);

    Matrix<Real, Eigen::Dynamic, 1> res;
    res.noalias() = interpolate(nodal_values, shapes);

    return res;
  }

  /// interpolate a field given the shape functions on the interpolation point
  template <class Derived1, class Derived2>
  static inline auto
  interpolate(const Eigen::MatrixBase<Derived1> &nodal_values,
              const Eigen::MatrixBase<Derived2> &shapes);

  /// interpolate a field given the shape functions on the interpolations points
  template <class Derived1, class Derived2, class Derived3>
  static inline void
  interpolate(const Eigen::MatrixBase<Derived1> &nodal_values,
              const Eigen::MatrixBase<Derived2> &Ns,
              const Eigen::MatrixBase<Derived3> &interpolated);

  /// compute the gradient of a given field on the given natural coordinates
  template <class D1, class D2, class D3>
  static inline void
  gradientOnNaturalCoordinates(const Eigen::MatrixBase<D1> &natural_coords,
                               const Eigen::MatrixBase<D2> &f,
                               const Eigen::MatrixBase<D3> &dfds);

public:
  static constexpr auto getShapeSize() {
    return InterpolationProperty<interpolation_type>::nb_nodes_per_element;
  }
  static constexpr auto getShapeDerivativesSize() {
    return (InterpolationProperty<interpolation_type>::nb_nodes_per_element *
            InterpolationProperty<interpolation_type>::natural_space_dimension);
  }
  static constexpr auto getNaturalSpaceDimension() {
    return InterpolationProperty<interpolation_type>::natural_space_dimension;
  }
  static constexpr auto getNbNodesPerInterpolationElement() {
    return InterpolationProperty<interpolation_type>::nb_nodes_per_element;
  }
};

/* -------------------------------------------------------------------------- */
/* Integration                                                                */
/* -------------------------------------------------------------------------- */
template <GaussIntegrationType git_class, Int nb_points>
struct GaussIntegrationTypeData {
  /// quadrature points in natural coordinates
  static Real quad_positions[];
  /// weights for the Gauss integration
  static Real quad_weights[];
};

template <ElementType type,
          Int n = ElementClassProperty<type>::polynomial_degree>
class GaussIntegrationElement {
  static constexpr InterpolationType itp_type =
      ElementClassProperty<type>::interpolation_type;
  using interpolation_property = InterpolationProperty<itp_type>;

public:
  static constexpr Int getNbQuadraturePoints();
  static constexpr auto getQuadraturePoints()
      -> Matrix<Real, interpolation_property::natural_space_dimension,
                GaussIntegrationElement::getNbQuadraturePoints()>;
  static constexpr auto getWeights()
      -> Vector<Real, GaussIntegrationElement::getNbQuadraturePoints()>;
};

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */
template <ElementType element_type,
          ElementKind element_kind =
              ElementClassProperty<element_type>::element_kind>
class ElementClass
    : public GeometricalElement<
          ElementClassProperty<element_type>::geometrical_type>,
      public InterpolationElement<
          ElementClassProperty<element_type>::interpolation_type> {
protected:
  using geometrical_element =
      GeometricalElement<ElementClassProperty<element_type>::geometrical_type>;
  using interpolation_element = InterpolationElement<
      ElementClassProperty<element_type>::interpolation_type>;

  using element_property = ElementClassProperty<element_type>;
  using interpolation_property =
      typename interpolation_element::interpolation_property;

public:
  /**
   * compute @f$ J = \frac{\partial x_j}{\partial s_i} @f$ the variation of real
   * coordinates along with variation of natural coordinates on a given point in
   * natural coordinates
   */
  template <class D1, class D2>
  static inline decltype(auto)
  computeJMat(const Eigen::MatrixBase<D1> &dnds,
              const Eigen::MatrixBase<D2> &node_coords);

  /**
   * compute the Jacobian matrix by computing the variation of real coordinates
   * along with variation of natural coordinates on a given set of points in
   * natural coordinates
   */
  template <class D>
  static inline void computeJMat(const Tensor3Base<Real> &dnds,
                                 const Eigen::MatrixBase<D> &node_coords,
                                 Tensor3Base<Real> &J);

  /// compute the jacobians of a serie of natural coordinates
  template <class D1, class D2, class D3>
  static inline void
  computeJacobian(const Eigen::MatrixBase<D1> &natural_coords,
                  const Eigen::MatrixBase<D2> &node_coords,
                  Eigen::MatrixBase<D3> &jacobians);

  /// compute jacobian (or integration variable change factor) for a set of
  /// points
  template <class D>
  static inline void computeJacobian(const Tensor3Base<Real> &J,
                                     Eigen::MatrixBase<D> &jacobians);

  /// compute jacobian (or integration variable change factor) for a given point
  template <class D>
  static inline Real computeJacobian(const Eigen::MatrixBase<D> &J);

  /// compute shape derivatives (input is dxds) for a set of points
  static inline void computeShapeDerivatives(const Tensor3Base<Real> &J,
                                             const Tensor3Base<Real> &dnds,
                                             Tensor3Base<Real> &shape_deriv);

  /// compute shape derivatives (input is dxds) for a given point
  template <class D1, class D2, class D3>
  static inline void
  computeShapeDerivatives(const Eigen::MatrixBase<D1> &J,
                          const Eigen::MatrixBase<D2> &dnds,
                          Eigen::MatrixBase<D3> &shape_deriv);

  /// compute the normal of a surface defined by the function f
  template <class D1, class D2, class D3>
  static inline void
  computeNormalsOnNaturalCoordinates(const Eigen::MatrixBase<D1> &coord,
                                     const Eigen::MatrixBase<D2> &f,
                                     Eigen::MatrixBase<D3> &normals);

  /// get natural coordinates from real coordinates
  template <class D1, class D2, class D3,
            aka::enable_if_vectors_t<D1, D3> * = nullptr>
  static inline void inverseMap(const Eigen::MatrixBase<D1> &real_coords,
                                const Eigen::MatrixBase<D2> &node_coords,
                                const Eigen::MatrixBase<D3> &natural_coords,
                                Int max_iterations = 100,
                                Real tolerance = 1e-10);

  /// get natural coordinates from real coordinates
  template <class D1, class D2, class D3,
            aka::enable_if_matrices_t<D1, D3> * = nullptr>
  static inline void inverseMap(const Eigen::MatrixBase<D1> &real_coords,
                                const Eigen::MatrixBase<D2> &node_coords,
                                const Eigen::MatrixBase<D3> &natural_coords_,
                                Int max_iterations = 100,
                                Real tolerance = 1e-10);

public:
  static constexpr auto getKind() { return element_kind; }
  static constexpr auto getSpatialDimension() {
    return ElementClassProperty<element_type>::spatial_dimension;
  }

  using element_class_extra_geom_property =
      ElementClassExtraGeometryProperties<element_type>;

  static constexpr decltype(auto) getP1ElementType() {
    return element_class_extra_geom_property::p1_type;
  }
  static constexpr decltype(auto) getFacetType(UInt t = 0) {
    return element_class_extra_geom_property::facet_type[t];
  }
  static constexpr decltype(auto) getFacetTypes();
};

/* -------------------------------------------------------------------------- */
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "geometrical_element_property.hh"
#include "interpolation_element_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include "element_class_tmpl.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_not_defined, _gt_not_defined,
                                     _itp_not_defined, _ek_not_defined, 0,
                                     _git_not_defined, 0);
} // namespace akantu

#include "element_class_hexahedron_8_inline_impl.hh"
#include "element_class_pentahedron_6_inline_impl.hh"
/* keep order */
#include "element_class_hexahedron_20_inline_impl.hh"
#include "element_class_pentahedron_15_inline_impl.hh"
#include "element_class_point_1_inline_impl.hh"
#include "element_class_quadrangle_4_inline_impl.hh"
#include "element_class_quadrangle_8_inline_impl.hh"
#include "element_class_segment_2_inline_impl.hh"
#include "element_class_segment_3_inline_impl.hh"
#include "element_class_tetrahedron_10_inline_impl.hh"
#include "element_class_tetrahedron_4_inline_impl.hh"
#include "element_class_triangle_3_inline_impl.hh"
#include "element_class_triangle_6_inline_impl.hh"

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "element_class_structural.hh"
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element.hh"
#endif

#if defined(AKANTU_IGFEM)
#include "element_class_igfem.hh"
#endif

#endif /* AKANTU_ELEMENT_CLASS_HH_ */
