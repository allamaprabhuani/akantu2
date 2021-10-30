/**
 * @file   element_class_bernoulli_beam_inline_impl.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Fri Feb 05 2021
 *
 * @brief  Specialization of the element_class class for the type
 * _bernoulli_beam_2
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

/**
 * @verbatim
   --x-----q1----|----q2-----x---> x
    -1          0            1
 @endverbatim
 *
 */

/* -------------------------------------------------------------------------- */
//#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_HH_
#define AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_2,
                                                     _itp_lagrange_segment_2, 3,
                                                     2, 6);

AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_bernoulli_beam_3,
                                                     _itp_lagrange_segment_2, 6,
                                                     4, 6);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_2,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_2,
                                                _segment_2, _ek_structural, 2,
                                                _git_segment, 3);

AKANTU_DEFINE_STRUCTURAL_ELEMENT_CLASS_PROPERTY(_bernoulli_beam_3,
                                                _gt_segment_2,
                                                _itp_bernoulli_beam_3,
                                                _segment_2, _ek_structural, 3,
                                                _git_segment, 3);

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeShapes(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & real_coord,
    Eigen::MatrixBase<D3> & N) {
  Eigen::Matrix<Real, 2, 1> L;
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeShapes(
      natural_coords, L);
  Eigen::Matrix<Real, 2, 4> H;
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
      natural_coords, real_coord, H);

  // clang-format off
  //    u1   v1      t1      u2   v2      t2
  N << L(0), 0      , 0      , L(1), 0      , 0      ,  // u
       0   , H(0, 0), H(0, 1), 0   , H(0, 2), H(0, 3),  // v
       0   , H(1, 0), H(1, 1), 0   , H(1, 2), H(1, 3); // theta
  // clang-format on
}

template <>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeShapes(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & real_coord,
    Eigen::MatrixBase<D3> & N) {
  Eigen::Matrix<Real, 2, 1> L;
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeShapes(
      natural_coords, L);
  Eigen::Matrix<Real, 2, 4> H;
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
      natural_coords, real_coord, H);

  // clang-format off
  //    u1    v1       w1       tx1   ty1       tz1    u2      v2       w2       tx2   ty2       tz2
  N << L(0), 0      , 0      , 0   , 0       , 0      , L(1), 0      , 0      , 0   , 0       , 0      ,  // u
       0   , H(0, 0), 0      , 0   , 0       , H(0, 1), 0   , H(0, 2), 0      , 0   , 0       , H(0, 3),  // v
       0   , 0      , H(0, 0), 0   , -H(0, 1), 0      , 0   , 0      , H(0, 2), 0   , -H(0, 3), 0      ,  // w
       0   , 0      , 0      , L(0), 0       , 0      , 0   , 0      , 0      , L(1), 0       , 0      ,  // thetax
       0   , 0      , H(1, 0), 0   , -H(1, 1), 0      , 0   , 0      , H(1, 2), 0   , -H(1, 3), 0      ,  // thetay
       0   , H(1, 0), 0      , 0   , 0       , H(1, 1), 0   , H(1, 2), 0      , 0   , 0       , H(1, 3); // thetaz
  // clang-format on
}

/* -------------------------------------------------------------------------- */
#if 0
template <>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeShapesDisplacements(
    const Vector<Real> & natural_coords, const Matrix<Real> & real_coord,
    Matrix<Real> & N) {
}
#endif

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2, class D3>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeDNDS(
    const Eigen::MatrixBase<D1> & Xs, const Eigen::MatrixBase<D2> & xs,
    Eigen::MatrixBase<D3> & dnds) {
  Eigen::Matrix<Real, 1, 2> L;
  InterpolationElement<_itp_lagrange_segment_2, _itk_lagrangian>::computeDNDS(
      Xs, L);
  Eigen::Matrix<Real, 1, 4> H;
  InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(Xs, xs, H);

  // Storing the derivatives in dnds
  dnds.block(0, 0, L.rows(), L.cols()) = L;
  dnds.block(0, 2, H.rows(), H.cols()) = H;
}

/* -------------------------------------------------------------------------- */
template <>
template<class D1, class D2>
inline void
InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::arrangeInVoigt(
    const Eigen::MatrixBase<D1> & dnds, Eigen::MatrixBase<D2> & B) {
  auto L = dnds.block(0, 0, 1, 2); // Lagrange shape derivatives
  auto H = dnds.block(0, 2, 1, 4); // Hermite shape derivatives
  // clang-format off
  //    u1       v1       t1        u2        v2        t2
  B << L(0, 0), 0,       0,        L(0, 1),  0,        0      ,
       0,      -H(0, 0), -H(0, 1), 0,       -H(0, 2), -H(0, 3);
  // clang-format on
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2, class D3>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::computeDNDS(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & real_coord, Eigen::MatrixBase <D3> &dnds) {
  InterpolationElement<_itp_bernoulli_beam_2, _itk_structural>::computeDNDS(
      natural_coords, real_coord, dnds);
}

/* -------------------------------------------------------------------------- */
template <>
template<class D1, class D2>
inline void
InterpolationElement<_itp_bernoulli_beam_3, _itk_structural>::arrangeInVoigt(
    const Eigen::MatrixBase<D1> & dnds, Eigen::MatrixBase<D2> & B) {
  auto L = dnds.block(0, 0, 1, 2); // Lagrange shape derivatives
  auto H = dnds.block(0, 2, 1, 4); // Hermite shape derivatives

  // clang-format off
  //    u1       v1        w1        x1       y1        z1        u2       v2        w2         x2       y2       z2
  B << L(0, 0), 0       , 0       , 0      , 0       , 0       , L(0, 1), 0       , 0        , 0      , 0        , 0      ,  // eps
       0      , -H(0, 0), 0       , 0      , 0       , -H(0, 1), 0      , -H(0, 2), 0        , 0      , 0        ,-H(0, 3),  // chi strong axis
       0      , 0       , -H(0, 0), 0      , H(0, 1) , 0       , 0      , 0       , -H(0, 2) , 0      , H(0, 3)  , 0      ,  // chi weak axis
       0      , 0       , 0       , L(0, 0), 0       , 0       , 0      , 0       , 0        , L(0, 1), 0        , 0      ; // chi torsion
  // clang-format on
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2, class D3>
inline void ElementClass<_bernoulli_beam_2>::computeRotationMatrix(
    Eigen::MatrixBase<D1> & R, const Eigen::MatrixBase<D2> & X,
    const Eigen::MatrixBase<D3> &) {
  auto && x2 = X(1); // X2
  auto && x1 = X(0); // X1

  auto cs = (x2 - x1) / (x2 - x1).norm();
  
  auto c = cs(0);
  auto s = cs(1);

  // clang-format off
  /// Definition of the rotation matrix
  R << c,  s,  0.,
      -s,  c,  0.,
       0., 0., 1.;
  // clang-format on
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2, class D3>
inline void ElementClass<_bernoulli_beam_3>::computeRotationMatrix(
    Eigen::MatrixBase<D1> & R, const Eigen::MatrixBase<D2> & X,
    const Eigen::MatrixBase<D3> & n) {
  Vector<Real> x2 = X(1); // X2
  Vector<Real> x1 = X(0); // X1
  auto dim = X.rows();
  Eigen::Matrix<Real, 1, 3> x = (x2 - x1), nv = n;

  x.normalize();
  auto x_n = x.cross(nv);

  Matrix<Real> Pe(dim, dim);
  Pe << 1., 0., 0., 0., -1., 0., 0., 0., 1.;

  Matrix<Real> Pg(dim, dim);
  Pg(0) = x;
  Pg(1) = x_n;
  Pg(2) = n;

  Pe *= Pg.inverse();

  R.zero();
  /// Definition of the rotation matrix
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      R(i + dim, j + dim) = R(i, j) = Pe(i, j);
    }
  }
}

} // namespace akantu
#endif /* AKANTU_ELEMENT_CLASS_BERNOULLI_BEAM_INLINE_IMPL_HH_ */
