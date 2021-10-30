/**
 * @file   element_class_kirchhoff_shell_inline_impl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 04 2014
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Element class Kirchhoff Shell
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
//#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH_
#define AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(
    _itp_discrete_kirchhoff_triangle_18, _itp_lagrange_triangle_3, 6, 6, 21);
AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(
    _discrete_kirchhoff_triangle_18, _gt_triangle_3,
    _itp_discrete_kirchhoff_triangle_18, _triangle_3, _ek_structural, 3,
    _git_triangle, 2);

/* -------------------------------------------------------------------------- */
namespace detail {
  template <class D>
  inline decltype(auto)
  computeBasisChangeMatrix(const Eigen::MatrixBase<D> & X) {
    auto && X1 = X(0);
    auto && X2 = X(1);
    auto && X3 = X(2);

    Eigen::Matrix<Real, 1, 3> a1 = X2 - X1;
    Eigen::Matrix<Real, 1, 3> a2 = X3 - X1;

    a1.normalized();
    auto && e3 = a1.cross(a2).normalized();
    auto && e2 = e3.cross(a1);

    Eigen::Matrix<Real, 3, 3> P;
    P(0) = a1;
    P(1) = e2;
    P(2) = e3;

    return P.transpose();
  }
} // namespace detail

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2, class D3>
inline void
ElementClass<_discrete_kirchhoff_triangle_18>::computeRotationMatrix(
    Eigen::MatrixBase<D1> & R, const Eigen::MatrixBase<D2> & X,
    const Eigen::MatrixBase<D3> &) {
  auto dim = X.rows();

  auto && P = detail::computeBasisChangeMatrix(X);

  R.zero();
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      R(i + dim, j + dim) = R(i, j) = P(i, j);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeShapes(
    const Eigen::MatrixBase<D1> & /*natural_coords*/,
    const Eigen::MatrixBase<D2> & /*real_coord*/,
    Eigen::MatrixBase<D3> & /*N*/) {}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18>::computeDNDS(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & real_coordinates, Eigen::MatrixBase<D3> & B) {

  auto && P = detail::computeBasisChangeMatrix(real_coordinates);

  auto && X = P * real_coordinates;
  auto && X1 = X(0);
  auto && X2 = X(1);
  auto && X3 = X(2);

  std::array<Vector<Real>, 3> A = {X2 - X1, X3 - X2, X1 - X3};
  std::array<Real, 3> L;
  std::array<Real, 3> C;
  std::array<Real, 3> S;

  // Setting all last coordinates to 0
  std::for_each(A.begin(), A.end(), [](auto & a) { a(2) = 0; });
  // Computing lengths
  std::transform(A.begin(), A.end(), L.begin(),
                 [](auto & a) { return a.norm(); });
  // Computing cosines
  std::transform(A.begin(), A.end(), L.begin(), C.begin(),
                 [](auto & a, auto & l) { return a(0) / l; });
  // Computing sines
  std::transform(A.begin(), A.end(), L.begin(), S.begin(),
                 [](auto & a, auto & l) { return a(1) / l; });

  // Natural coordinates
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);

  // Derivative of quadratic interpolation functions
  Eigen::Matrix<Real, 2, 3> dP;
  dP << 4 * (1 - 2 * xi - eta), 4 * eta, -4 * eta, -4 * xi, 4 * xi,
      4 * (1 - xi - 2 * eta);

  Eigen::Matrix<Real, 2, 3> dNx1;
  dNx1 << 3. / 2 * (dP(0, 0) * C[0] / L[0] - dP(0, 2) * C[2] / L[2]),
      3. / 2 * (dP(0, 1) * C[1] / L[1] - dP(0, 0) * C[0] / L[0]),
      3. / 2 * (dP(0, 2) * C[2] / L[2] - dP(0, 1) * C[1] / L[1]),
      3. / 2 * (dP(1, 0) * C[0] / L[0] - dP(1, 2) * C[2] / L[2]),
      3. / 2 * (dP(1, 1) * C[1] / L[1] - dP(1, 0) * C[0] / L[0]),
      3. / 2 * (dP(1, 2) * C[2] / L[2] - dP(1, 1) * C[1] / L[1]);
  Eigen::Matrix<Real, 2, 3> dNx2;
  dNx2 << -1 - 3. / 4 * (dP(0, 0) * C[0] * C[0] + dP(0, 2) * C[2] * C[2]),
      1 - 3. / 4 * (dP(0, 1) * C[1] * C[1] + dP(0, 0) * C[0] * C[0]),
      -3. / 4 * (dP(0, 2) * C[2] * C[2] + dP(0, 1) * C[1] * C[1]),
      -1 - 3. / 4 * (dP(1, 0) * C[0] * C[0] + dP(1, 2) * C[2] * C[2]),
      -3. / 4 * (dP(1, 1) * C[1] * C[1] + dP(1, 0) * C[0] * C[0]),
      1 - 3. / 4 * (dP(1, 2) * C[2] * C[2] + dP(1, 1) * C[1] * C[1]);

  Eigen::Matrix<Real, 2, 3> dNx3;
  dNx3 << -3. / 4 * (dP(0, 0) * C[0] * S[0] + dP(0, 2) * C[2] * S[2]),
      -3. / 4 * (dP(0, 1) * C[1] * S[1] + dP(0, 0) * C[0] * S[0]),
      -3. / 4 * (dP(0, 2) * C[2] * S[2] + dP(0, 1) * C[1] * S[1]),
      -3. / 4 * (dP(1, 0) * C[0] * S[0] + dP(1, 2) * C[2] * S[2]),
      -3. / 4 * (dP(1, 1) * C[1] * S[1] + dP(1, 0) * C[0] * S[0]),
      -3. / 4 * (dP(1, 2) * C[2] * S[2] + dP(1, 1) * C[1] * S[1]);
  Eigen::Matrix<Real, 2, 3> dNy1;
  dNy1 << 3. / 2 * (dP(0, 0) * S[0] / L[0] - dP(0, 2) * S[2] / L[2]),
      3. / 2 * (dP(0, 1) * S[1] / L[1] - dP(0, 0) * S[0] / L[0]),
      3. / 2 * (dP(0, 2) * S[2] / L[2] - dP(0, 1) * S[1] / L[1]),
      3. / 2 * (dP(1, 0) * S[0] / L[0] - dP(1, 2) * S[2] / L[2]),
      3. / 2 * (dP(1, 1) * S[1] / L[1] - dP(1, 0) * S[0] / L[0]),
      3. / 2 * (dP(1, 2) * S[2] / L[2] - dP(1, 1) * S[1] / L[1]);
  auto dNy2 = dNx3;

  Eigen::Matrix<Real, 2, 3> dNy3;
  dNy3 << -1 - 3. / 4 * (dP(0, 0) * S[0] * S[0] + dP(0, 2) * S[2] * S[2]),
      1 - 3. / 4 * (dP(0, 1) * S[1] * S[1] + dP(0, 0) * S[0] * S[0]),
      -3. / 4 * (dP(0, 2) * S[2] * S[2] + dP(0, 1) * S[1] * S[1]),
      -1 - 3. / 4 * (dP(1, 0) * S[0] * S[0] + dP(1, 2) * S[2] * S[2]),
      -3. / 4 * (dP(1, 1) * S[1] * S[1] + dP(1, 0) * S[0] * S[0]),
      1 - 3. / 4 * (dP(1, 2) * S[2] * S[2] + dP(1, 1) * S[1] * S[1]);

  // Derivative of linear (membrane mode) functions
  Eigen::Matrix<Real, 2, 3> dNm;
  InterpolationElement<_itp_lagrange_triangle_3, _itk_lagrangian>::computeDNDS(
      natural_coords, dNm);

  UInt i = 0;
  for (auto && mat : {dNm, dNx1, dNx2, dNx3, dNy1, dNy2, dNy3}) {
    B.block(2, 3, 0, i) = mat;
    i += mat.cols();
  }
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2>
inline void
InterpolationElement<_itp_discrete_kirchhoff_triangle_18, _itk_structural>::
    arrangeInVoigt(const Eigen::MatrixBase<D1> & dnds,
                   Eigen::MatrixBase<D2> & B) {
  Eigen::Matrix<Real, 2, 3> dNm, dNx1, dNx2, dNx3, dNy1, dNy2, dNy3;

  UInt i = 0;
  for (auto && mat : {&dNm, &dNx1, &dNx2, &dNx3, &dNy1, &dNy2, &dNy3}) {
    *mat = dnds.block(2, 3, 0, i);
    i += mat->cols();
  }

  for (UInt i = 0; i < 3; ++i) {
    // clang-format off
    Eigen::Matrix<Real, 3, 6> Bm;
    Bm << dNm(0, i), 0,         0, 0, 0, 0,
          0,         dNm(1, i), 0, 0, 0, 0,
          dNm(1, i), dNm(0, i), 0, 0, 0, 0;
    Eigen::Matrix<Real, 3, 6> Bf;
    Bf << 0, 0, dNx1(0, i),              -dNx3(0, i),              dNx2(0, i),              0,
          0, 0, dNy1(1, i),              -dNy3(1, i),              dNy2(1, i),              0,
          0, 0, dNx1(1, i) + dNy1(0, i), -dNx3(1, i) - dNy3(0, i), dNx2(1, i) + dNy2(0, i), 0;
    // clang-format on
    B.block(3, 6, 0, i * 6) = Bm;
    B.block(3, 6, 3, i * 6) = Bf;
  }
}

} // namespace akantu

#endif /* AKANTU_ELEMENT_CLASS_KIRCHHOFF_SHELL_INLINE_IMPL_HH_ */
