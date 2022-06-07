/**
 * @file   element_class_structural.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Mon Feb 01 2021
 *
 * @brief  Specialization of the element classes for structural elements
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
#include "aka_iterators.hh"
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_
#define AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_

namespace akantu {

/// Macro to generate the InterpolationProperty structures for different
/// interpolation types
#define AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(                  \
    itp_type, itp_geom_type, ndof, nb_stress, nb_dnds_cols)                    \
  template <> struct InterpolationProperty<itp_type> {                         \
    static constexpr InterpolationKind kind{_itk_structural};                  \
    static constexpr Int nb_nodes_per_element{                                 \
        InterpolationProperty<itp_geom_type>::nb_nodes_per_element};           \
    static constexpr InterpolationType itp_geometry_type{itp_geom_type};       \
    static constexpr Int natural_space_dimension{                              \
        InterpolationProperty<itp_geom_type>::natural_space_dimension};        \
    static constexpr Int nb_degree_of_freedom{ndof};                           \
    static constexpr Int nb_stress_components{nb_stress};                      \
    static constexpr Int dnds_columns{nb_dnds_cols};                           \
  }

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type>
class InterpolationElement<interpolation_type, _itk_structural> {
public:
  using interpolation_property = InterpolationProperty<interpolation_type>;

  /// compute the shape values for a given point in natural coordinates
  template <class D1, class D2, class D3>
  static inline void computeShapes(const Eigen::MatrixBase<D1> & natural_coord,
                                   const Eigen::MatrixBase<D2> & real_coord,
                                   Eigen::MatrixBase<D3> & N);

  /// compute the shape values for a given set of points in natural coordinates
  template <class D1, class D2, class D3>
  static inline void computeShapes(const Eigen::MatrixBase<D1> & Xs,
                                   const Eigen::MatrixBase<D2> & x,
                                   const Eigen::MatrixBase<D3> & T,
                                   TensorBase<Real, 3> & Ns) {
    Matrix<Real> N(Ns.size(1), Ns.size(2));
    for (auto && data : zip(Xs, Ns)) {
      auto && X = std::get<0>(data);
      auto && N_T = std::get<1>(data);

      computeShapes(X, x, N);
      N_T = N.block(0, 0, N_T.rows(), N_T.cols()) * T;
    }
  }

  template <class D1, class D2, class D3>
  static inline void computeShapesMass(const Eigen::MatrixBase<D1> & Xs,
                                       const Eigen::MatrixBase<D2> & x,
                                       const Eigen::MatrixBase<D3> & T,
                                       TensorBase<Real, 3> & Ns) {
    for (int i = 0; i < Xs.cols(); ++i) {
      auto N_T = Ns(i);

      Matrix<Real> N(interpolation_property::nb_degree_of_freedom, N_T.cols());
      computeShapes(Xs(i), x, N);

      N_T = N.block(0, 0, N_T.rows(), N_T.cols()) * T;
    }
  }

  /// compute shape derivatives (input is dxds) for a set of points
  template <class D>
  static inline void computeShapeDerivatives(const TensorBase<Real, 3> & Js,
                                             const TensorBase<Real, 3> & DNDSs,
                                             const Eigen::MatrixBase<D> & R,
                                             TensorBase<Real, 3> & Bs) {
    for (Int i = 0; i < Js.size(2); ++i) {
      auto && DNDX = Js(i).inverse() * DNDSs(i);
      auto && B_R = Bs(i);
      Matrix<Real> B(B_R.rows(), B_R.cols());
      arrangeInVoigt(DNDX, B);
      B_R = B * R;
    }
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with variation of natural coordinates on a given set
   * of points in natural coordinates
   */
  template <typename D1, typename D2>
  static inline void computeDNDS(const Eigen::MatrixBase<D1> & Xs,
                                 const Eigen::MatrixBase<D2> & xs,
                                 TensorBase<Real, 3> & dnds) {
    for (auto && data : zip(Xs, dnds))
      computeDNDS(std::get<0>(data), xs, std::get<1>(data));
  }

  /**
   * compute @f$ B_{ij} = \frac{\partial N_j}{\partial S_i} @f$ the variation of
   * shape functions along with
   * variation of natural coordinates on a given point in natural
   * coordinates
   */
  template <typename D1, typename D2, typename D3>
  static inline void computeDNDS(const Eigen::MatrixBase<D1> & Xs,
                                 const Eigen::MatrixBase<D2> & xs,
                                 Eigen::MatrixBase<D3> & dnds);

  /**
   * arrange B in Voigt notation from DNDS
   */
  template <class D1, class D2>
  static inline void arrangeInVoigt(const Eigen::MatrixBase<D1> & dnds,
                                    Eigen::MatrixBase<D2> & B) {
    // Default implementation assumes dnds is already in Voigt notation
    B = dnds;
  }

public:
  static inline constexpr auto getNbNodesPerInterpolationElement() {
    return interpolation_property::nb_nodes_per_element;
  }

  static inline constexpr auto getShapeSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_degree_of_freedom;
  }
  static inline constexpr auto getShapeIndependantSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_stress_components;
  }
  static inline constexpr auto getShapeDerivativesSize() {
    return interpolation_property::nb_nodes_per_element *
           interpolation_property::nb_degree_of_freedom *
           interpolation_property::nb_stress_components;
  }
  static inline constexpr auto getNaturalSpaceDimension() {
    return interpolation_property::natural_space_dimension;
  }
  static inline constexpr auto getNbDegreeOfFreedom() {
    return interpolation_property::nb_degree_of_freedom;
  }
  static inline constexpr auto getNbStressComponents() {
    return interpolation_property::nb_stress_components;
  }
};

/// Macro to generate the element class structures for different structural
/// element types
/* -------------------------------------------------------------------------- */
#define AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(                       \
    elem_type, geom_type, interp_type, parent_el_type, sp, gauss_int_type,     \
    min_int_order)                                                             \
  template <> struct ElementClassProperty<elem_type> {                         \
    static constexpr GeometricalType geometrical_type{geom_type};              \
    static constexpr InterpolationType interpolation_type{interp_type};        \
    static constexpr ElementType parent_element_type{parent_el_type};          \
    static constexpr ElementKind element_kind{_ek_structural};                 \
    static constexpr Int spatial_dimension{sp};                                \
    static constexpr GaussIntegrationType gauss_integration_type{              \
        gauss_int_type};                                                       \
    static constexpr Int polynomial_degree{min_int_order};                     \
  }

/* -------------------------------------------------------------------------- */
/* ElementClass for structural elements                                       */
/* -------------------------------------------------------------------------- */
template <ElementType element_type>
class ElementClass<element_type, _ek_structural>
    : public GeometricalElement<
          ElementClassProperty<element_type>::geometrical_type>,
      public InterpolationElement<
          ElementClassProperty<element_type>::interpolation_type> {
protected:
  using geometrical_element =
      GeometricalElement<ElementClassProperty<element_type>::geometrical_type>;
  using interpolation_element = InterpolationElement<
      ElementClassProperty<element_type>::interpolation_type>;
  using parent_element =
      ElementClass<ElementClassProperty<element_type>::parent_element_type>;

public:
  template <class D1, class D2, class D3>
  static inline void
  computeRotationMatrix(Eigen::MatrixBase<D1> & /*R*/,
                        const Eigen::MatrixBase<D2> & /*X*/,
                        const Eigen::MatrixBase<D3> & /*extra_normal*/) {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute jacobian (or integration variable change factor) for a given point
  template <typename D1, typename D2, typename D3>
  static inline void computeJMat(const Eigen::MatrixBase<D1> & natural_coords,
                                 const Eigen::MatrixBase<D2> & Xs,
                                 Eigen::MatrixBase<D3> & J) {
    Matrix<Real> dnds(Xs.rows(), Xs.cols());
    parent_element::computeDNDS(natural_coords, dnds);
    J = dnds * Xs.transpose();
  }

  template <typename D1, typename D2>
  static inline void computeJMat(const Eigen::MatrixBase<D1> & Xs,
                                 const Eigen::MatrixBase<D2> & xs,
                                 Tensor3<Real> & Js) {
    for (auto && data : zip(Xs, Js)) {
      computeJMat(std::get<0>(data), xs, std::get<1>(data));
    }
  }

  template <typename D1, typename D2, typename D3,
            std::enable_if_t<aka::is_vector<D3>::value> * = nullptr>
  static inline void computeJacobian(const Eigen::MatrixBase<D1> & Xs,
                                     const Eigen::MatrixBase<D2> & xs,
                                     Eigen::MatrixBase<D3> & jacobians) {
    using itp = typename interpolation_element::interpolation_property;
    Tensor3<Real> Js(itp::natural_space_dimension, itp::natural_space_dimension,
                     Xs.cols());
    computeJMat(Xs, xs, Js);
    for (auto && data : zip(jacobians, Js)) {
      std::get<0>(data) = std::get<1>(data).determinant();
    }
  }

  template <typename D1, typename D2>
  static inline void computeRotation(const Eigen::MatrixBase<D1> & xs,
                                     Eigen::MatrixBase<D2> & R);

public:
  static constexpr AKANTU_GET_MACRO_AUTO_NOT_CONST(Kind, _ek_structural);
  static constexpr AKANTU_GET_MACRO_AUTO_NOT_CONST(P1ElementType, _not_defined);
  static constexpr AKANTU_GET_MACRO_AUTO_NOT_CONST(FacetType, _not_defined);
  static constexpr auto getFacetType(__attribute__((unused)) Int t = 0) {
    return _not_defined;
  }
  static constexpr AKANTU_GET_MACRO_AUTO_NOT_CONST(
      SpatialDimension, ElementClassProperty<element_type>::spatial_dimension);
  static constexpr auto getFacetTypes() {
    return ElementClass<_not_defined>::getFacetTypes();
  }
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "element_class_hermite_inline_impl.hh"
/* keep order */
#include "element_class_bernoulli_beam_inline_impl.hh"
#include "element_class_kirchhoff_shell_inline_impl.hh"
/* -------------------------------------------------------------------------- */

#endif /* AKANTU_ELEMENT_CLASS_STRUCTURAL_HH_ */
