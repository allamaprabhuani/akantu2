/**
 * @file   element_class_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Implementation of the inline templated function of the element class
 * descriptions
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
//#include "element_class.hh"
#include "aka_iterators.hh"
#include "gauss_integration_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_TMPL_HH_
#define AKANTU_ELEMENT_CLASS_TMPL_HH_

namespace akantu {

template <ElementType element_type, ElementKind element_kind>
inline constexpr decltype(auto)
ElementClass<element_type, element_kind>::getFacetTypes() {
  return Eigen::Map<const Eigen::Matrix<
      ElementType, geometrical_element::getNbFacetTypes(), 1>>(
      element_class_extra_geom_property::facet_type.data());
}

/* -------------------------------------------------------------------------- */
/* GeometricalElement                                                         */
/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
template <Idx t, std::size_t size, std::enable_if_t<not(t < size)> *>
inline constexpr decltype(auto)
GeometricalElement<geometrical_type,
                   shape>::getFacetLocalConnectivityPerElement() {
  AKANTU_ERROR("Not a valid facet id for this element type");
}

template <GeometricalType geometrical_type, GeometricalShapeType shape>
template <Idx t, std::size_t size, std::enable_if_t<(t < size)> *>
inline constexpr decltype(auto)
GeometricalElement<geometrical_type,
                   shape>::getFacetLocalConnectivityPerElement() {
  Int pos = 0;
  for (Int i = 0; i < t; ++i) {
    pos += geometrical_property::nb_facets[i] *
           geometrical_property::nb_nodes_per_facet[i];
  }

  return Eigen::Map<
      const Eigen::Matrix<Idx, geometrical_property::nb_facets[t],
                          geometrical_property::nb_nodes_per_facet[t]>>(
      geometrical_property::facet_connectivity_vect.data() + pos);
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline constexpr decltype(auto)
GeometricalElement<geometrical_type,
                   shape>::getFacetLocalConnectivityPerElement(Idx t) {

  Int pos = 0;
  for (Int i = 0; i < t; ++i) {
    pos += geometrical_property::nb_facets[i] *
           geometrical_property::nb_nodes_per_facet[i];
  }

  return Eigen::Map<const Eigen::Matrix<Idx, Eigen::Dynamic, Eigen::Dynamic>>(
      geometrical_property::facet_connectivity_vect.data() + pos,
      geometrical_property::nb_facets[t],
      geometrical_property::nb_nodes_per_facet[t]);
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline constexpr Int
GeometricalElement<geometrical_type, shape>::getNbFacetsPerElement() {
  Int total_nb_facets = 0;
  for (Int n = 0; n < geometrical_property::nb_facet_types; ++n) {
    total_nb_facets += geometrical_property::nb_facets[n];
  }

  return total_nb_facets;
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline constexpr Int
GeometricalElement<geometrical_type, shape>::getNbFacetsPerElement(Idx t) {
  return geometrical_property::nb_facets[t];
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
template <class D>
inline bool GeometricalElement<geometrical_type, shape>::contains(
    const Eigen::MatrixBase<D> & coords) {
  return GeometricalShapeContains<shape>::contains(coords);
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline bool GeometricalShapeContains<_gst_point>::contains(
    const Eigen::MatrixBase<D> & coords) {
  return (coords(0) < std::numeric_limits<Real>::epsilon());
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline bool GeometricalShapeContains<_gst_square>::contains(
    const Eigen::MatrixBase<D> & coords) {
  bool in = true;
  for (UInt i = 0; i < coords.size() && in; ++i) {
    in &= ((coords(i) >= -(1. + std::numeric_limits<Real>::epsilon())) &&
           (coords(i) <= (1. + std::numeric_limits<Real>::epsilon())));
  }
  return in;
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline bool GeometricalShapeContains<_gst_triangle>::contains(
    const Eigen::MatrixBase<D> & coords) {
  bool in = true;
  Real sum = 0;
  for (Int i = 0; (i < coords.size()) && in; ++i) {
    in &= ((coords(i) >= -(Math::getTolerance())) &&
           (coords(i) <= (1. + Math::getTolerance())));
    sum += coords(i);
  }
  if (in) {
    return (in && (sum <= (1. + Math::getTolerance())));
  }
  return in;
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline bool GeometricalShapeContains<_gst_prism>::contains(
    const Eigen::MatrixBase<D> & coords) {
  bool in = ((coords(0) >= -1.) && (coords(0) <= 1.)); // x in segment [-1, 1]

  // y and z in triangle
  in &= ((coords(1) >= 0) && (coords(1) <= 1.));
  in &= ((coords(2) >= 0) && (coords(2) <= 1.));
  Real sum = coords(1) + coords(2);

  return (in && (sum <= 1));
}

/* -------------------------------------------------------------------------- */
/* InterpolationElement                                                       */
/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type, InterpolationKind kind>
template <typename D1, typename D2,
          aka::enable_if_t<aka::are_matrices<D1, D2>::value> *>
inline void InterpolationElement<interpolation_type, kind>::computeShapes(
    const Eigen::MatrixBase<D1> & Xs, const Eigen::MatrixBase<D2> & N_) {

  Eigen::MatrixBase<D2> & N = const_cast<Eigen::MatrixBase<D2> &>(
      N_); // as advised by the Eigen developers

  for (auto && data : zip(Xs, N)) {
    computeShapes(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type, InterpolationKind kind>
template <class D>
inline void InterpolationElement<interpolation_type, kind>::computeDNDS(
    const Eigen::MatrixBase<D> & Xs, Tensor3Base<Real> & dNdS) {
  for (auto && data : zip(Xs, dNdS)) {
    computeDNDS(std::get<0>(data), std::get<1>(data));
  }
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate on a point a field for which values are given on the
 * node of the element using the shape functions at this interpolation point
 *
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param shapes value of shape functions at the interpolation point
 * @param interpolated interpolated value of f @f$ f_j(\xi) = \sum_i f_{n_i
 j} *N_i @f$
 */
template <InterpolationType interpolation_type, InterpolationKind kind>
template <typename Derived1, typename Derived2>
inline auto InterpolationElement<interpolation_type, kind>::interpolate(
    const Eigen::MatrixBase<Derived1> & nodal_values,
    const Eigen::MatrixBase<Derived2> & shapes) {
  return nodal_values * shapes;
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate on several points a field  for which values are given on the
 * node of the element using the shape functions at the interpolation point
 *
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param shapes value of shape functions at the interpolation point
 * @param interpolated interpolated values of f @f$ f_j(\xi) = \sum_i f_{n_i j}
 *N_i @f$
 */
template <InterpolationType interpolation_type, InterpolationKind kind>
template <typename Derived1, typename Derived2, typename Derived3>
inline void InterpolationElement<interpolation_type, kind>::interpolate(
    const Eigen::MatrixBase<Derived1> & nodal_values,
    const Eigen::MatrixBase<Derived2> & Ns,
    const Eigen::MatrixBase<Derived3> & interpolated_) {

  auto && interpolated = const_cast<Eigen::MatrixBase<Derived3> &>(
      interpolated_); // as advised by the Eigen developers

  auto nb_points = Ns.cols();
  for (auto p = 0; p < nb_points; ++p) {
    interpolated.col(p) = interpolate(nodal_values, Ns.col(p));
  }
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate the field on a point given in natural coordinates the field which
 * values are given on the node of the element
 *
 * @param natural_coords natural coordinates of point where to interpolate \xi
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param interpolated interpolated value of f @f$ f_j(\xi) = \sum_i f_{n_i j}
 *N_i @f$
 */
// template <InterpolationType interpolation_type, InterpolationKind kind>
// inline decltype(auto)
// InterpolationElement<interpolation_type,
// kind>::interpolateOnNaturalCoordinates(
//     const Ref<const VectorXr> & natural_coords,
//     const Ref<const MatrixXr> & nodal_values, Ref<VectorXr> interpolated) {
//   using interpolation = InterpolationProperty<interpolation_type>;
//   Eigen::Matrix<Real, interpolation::nb_nodes_per_element, 1> shapes;
//   computeShapes(natural_coords, shapes);

//   return interpolate(nodal_values, shapes);
// }

/* -------------------------------------------------------------------------- */
/// @f$ gradient_{ij} = \frac{\partial f_j}{\partial s_i} = \sum_k
/// \frac{\partial N_k}{\partial s_i}f_{j n_k} @f$
template <InterpolationType interpolation_type, InterpolationKind kind>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<interpolation_type, kind>::gradientOnNaturalCoordinates(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & f, const Eigen::MatrixBase<D3> & dfds_) {

  constexpr auto nsp =
      InterpolationProperty<interpolation_type>::natural_space_dimension;
  constexpr auto nnodes =
      InterpolationProperty<interpolation_type>::nb_nodes_per_element;
  Eigen::Matrix<Real, D3::ColsAtCompileTime, D2::ColsAtCompileTime> dnds(
      nsp, nnodes);
  auto & dfds = const_cast<Eigen::MatrixBase<D3> &>(dfds_);
  computeDNDS(natural_coords, dnds);
  dfds = f * dnds.transpose();
}

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D1, class D2>
inline decltype(auto) ElementClass<type, kind>::computeJMat(
    const Eigen::MatrixBase<D1> & dnds,
    const Eigen::MatrixBase<D2> & node_coords) {
  /// @f$ J = dxds = dnds * x @f$
  return dnds * node_coords.transpose();
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D>
inline void
ElementClass<type, kind>::computeJMat(const Tensor3Base<Real> & dnds,
                                      const Eigen::MatrixBase<D> & node_coords,
                                      Tensor3Base<Real> & J) {
  for (auto && data : zip(J, dnds)) {
    std::get<0>(data) = computeJMat(std::get<1>(data), node_coords);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D1, class D2, class D3>
inline void ElementClass<type, kind>::computeJacobian(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & node_coords,
    Eigen::MatrixBase<D3> & jacobians) {
  auto nb_points = natural_coords.cols();
  Matrix<Real, interpolation_property::natural_space_dimension,
         interpolation_property::nb_nodes_per_element>
      dnds;
  Matrix<Real> J(natural_coords.rows(), node_coords.rows());

  for (Int p = 0; p < nb_points; ++p) {
    interpolation_element::computeDNDS(natural_coords.col(p), dnds);
    J = computeJMat(dnds, node_coords);
    jacobians(p) = computeJacobian(J);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D>
inline void
ElementClass<type, kind>::computeJacobian(const Tensor3Base<Real> & J,
                                          Eigen::MatrixBase<D> & jacobians) {
  auto nb_points = J.size(2);
  for (Int p = 0; p < nb_points; ++p) {
    computeJacobian(J(p), jacobians.col(p));
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D>
inline Real
ElementClass<type, kind>::computeJacobian(const Eigen::MatrixBase<D> & J) {
  if (J.rows() == J.cols()) {
    return J.determinant();
  } else {
    switch (interpolation_property::natural_space_dimension) {
    case 1: {
      return J.norm();
    }
    case 2: {
      auto Jstatic =
          Eigen::Map<const Eigen::Matrix<Real, 2, 3>>(J.derived().data());
      return (Jstatic.row(0)).cross(Jstatic.row(1)).norm();
    }
    }
  }
  return 0; // avoids a warning
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D1, class D2, class D3>
inline void ElementClass<type, kind>::computeShapeDerivatives(
    const Eigen::MatrixBase<D1> & J, const Eigen::MatrixBase<D2> & dnds,
    Eigen::MatrixBase<D3> & shape_deriv) {
  shape_deriv = J.inverse() * dnds;
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeShapeDerivatives(
    const Tensor3Base<Real> & J, const Tensor3Base<Real> & dnds,
    Tensor3Base<Real> & shape_deriv) {
  auto nb_points = J.size(2);
  for (Int p = 0; p < nb_points; ++p) {
    auto && J_ = J(p);
    auto && dnds_ = dnds(p);
    auto && dndx_ = shape_deriv(p);
    dndx_ = J_.inverse() * dnds_;
    //    computeShapeDerivatives(J_, dnds_, dndx_);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D1, class D2, class D3>
inline void ElementClass<type, kind>::computeNormalsOnNaturalCoordinates(
    const Eigen::MatrixBase<D1> & coord, const Eigen::MatrixBase<D2> & f,
    Eigen::MatrixBase<D3> & normals) {
  auto dimension = normals.rows();
  auto nb_points = coord.cols();

  constexpr auto ndim = interpolation_property::natural_space_dimension;
  AKANTU_DEBUG_ASSERT((dimension - 1) == ndim,
                      "cannot extract a normal because of dimension mismatch "
                          << dimension - 1 << " " << ndim);

  Matrix<Real, Eigen::Dynamic, ndim> J(dimension, ndim);
  for (Int p = 0; p < nb_points; ++p) {
    interpolation_element::gradientOnNaturalCoordinates(coord.col(p), f, J);
    if (dimension == 2) {
      Math::normal2(J.data(), normals.col(p).data());
    }
    if (dimension == 3) {

      Math::normal3(J.col(0).data(), J.col(1).data(), normals.col(p).data());
    }
  }
}

/* ------------------------------------------------------------------------- */
/**
 * In the non linear cases we need to iterate to find the natural coordinates
 *@f$\xi@f$
 * provided real coordinates @f$x@f$.
 *
 * We want to solve: @f$ x- \phi(\xi) = 0@f$ with @f$\phi(\xi) = \sum_I N_I(\xi)
 *x_I@f$
 * the mapping function which uses the nodal coordinates @f$x_I@f$.
 *
 * To that end we use the Newton method and the following series:
 *
 * @f$ \frac{\partial \phi(x_k)}{\partial \xi} \left( \xi_{k+1} - \xi_k \right)
 *= x - \phi(x_k)@f$
 *
 * When we consider elements embedded in a dimension higher than them (2D
 *triangle in a 3D space for example)
 * @f$ J = \frac{\partial \phi(\xi_k)}{\partial \xi}@f$ is of dimension
 *@f$dim_{space} \times dim_{elem}@f$ which
 * is not invertible in most cases. Rather we can solve the problem:
 *
 * @f$ J^T J \left( \xi_{k+1} - \xi_k \right) = J^T \left( x - \phi(\xi_k)
 *\right) @f$
 *
 * So that
 *
 * @f$ d\xi = \xi_{k+1} - \xi_k = (J^T J)^{-1} J^T \left( x - \phi(\xi_k)
 *\right) @f$
 *
 * So that if the series converges we have:
 *
 * @f$ 0 = J^T \left( \phi(\xi_\infty) - x \right) @f$
 *
 * And we see that this is ill-posed only if @f$ J^T x = 0@f$ which means that
 *the vector provided
 * is normal to any tangent which means it is outside of the element itself.
 *
 * @param real_coords: the real coordinates the natural coordinates are sought
 *for
 * @param node_coords: the coordinates of the nodes forming the element
 * @param natural_coords: output->the sought natural coordinates
 * @param spatial_dimension: spatial dimension of the problem
 *
 **/
template <ElementType type, ElementKind kind>
template <class D1, class D2, class D3, aka::enable_if_vectors_t<D1, D3> *>
inline void ElementClass<type, kind>::inverseMap(
    const Eigen::MatrixBase<D1> & real_coords,
    const Eigen::MatrixBase<D2> & node_coords,
    const Eigen::MatrixBase<D3> & natural_coords_, Int max_iterations,
    Real tolerance) {
  auto & natural_coords = const_cast<Eigen::MatrixBase<D3> &>(
      natural_coords_); // as advised by the Eigen developers

  auto spatial_dimension = real_coords.size();
  auto dimension = natural_coords.size();

  // matrix copy of the real_coords
  //  MatrixProxy<const Real> mreal_coords(real_coords.data(),
  //  spatial_dimension,
  // 1);

  // initial guess
  natural_coords.zero();

  // real space coordinates provided by initial guess
  Matrix<Real> physical_guess(spatial_dimension, 1);

  // objective function f = real_coords - physical_guess
  Matrix<Real> f(spatial_dimension, 1);

  // J Jacobian matrix computed on the natural_guess
  Matrix<Real> J(dimension, spatial_dimension);

  // J^t
  Matrix<Real> Jt(spatial_dimension, dimension);

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  // do interpolation
  auto update_f = [&f, &physical_guess, &natural_coords, &node_coords,
                   &real_coords, spatial_dimension]() {
    Vector<Real> physical_guess_v(physical_guess.data(), spatial_dimension);
    interpolation_element::interpolateOnNaturalCoordinates(
        natural_coords, node_coords, physical_guess_v);

    // compute initial objective function value f = real_coords -
    // physical_guess
    f = real_coords - physical_guess;

    // compute initial error
    auto error = f.norm();
    return error;
  };

  auto inverse_map_error = update_f();
  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  UInt iterations{0};
  while (tolerance < inverse_map_error and iterations < max_iterations) {
    // compute J^t
    interpolation_element::gradientOnNaturalCoordinates(natural_coords,
                                                        node_coords, Jt);
    J = Jt.transpose();

    // compute G
    auto && G = J * J.transpose();

    // compute F
    auto && F = J.transpose() * G.inverse();

    // compute increment
    auto && dxi = F.transpose() * f;

    // update our guess
    natural_coords += dxi.col(0);

    inverse_map_error = update_f();
    iterations++;
  }

  if (iterations >= max_iterations) {
    AKANTU_EXCEPTION("The solver in inverse map did not converge");
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
template <class D1, class D2, class D3, aka::enable_if_matrices_t<D1, D3> *>
inline void ElementClass<type, kind>::inverseMap(
    const Eigen::MatrixBase<D1> & real_coords,
    const Eigen::MatrixBase<D2> & node_coords,
    const Eigen::MatrixBase<D3> & natural_coords_, Int max_iterations,
    Real tolerance) {
  auto & natural_coords = const_cast<Eigen::MatrixBase<D2> &>(
      natural_coords_); // as advised by the Eigen developers

  auto nb_points = real_coords.cols();
  for (Int p = 0; p < nb_points; ++p) {
    inverseMap(real_coords(p), node_coords, natural_coords(p), max_iterations,
               tolerance);
  }
}

} // namespace akantu

#endif /* AKANTU_ELEMENT_CLASS_TMPL_HH_ */
