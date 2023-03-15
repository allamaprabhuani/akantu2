/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/**
 * @verbatim
         \eta
      ^
 (-1,1)   |   (1,1)
     x---------x
     |    |    |
     |    |    |
   --|---------|----->  \xi
     |    |    |
     |    |    |
     x---------x
 (-1,-1)  |   (1,-1)
 @endverbatim
 *
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = - (1 - \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi) / 4 \\
 * N2 = (1 + \xi) (1 - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi) / 4 \\
 * N3 = (1 + \xi) (1 + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi) / 4 \\
 * N4 = (1 - \xi) (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = - (1 + \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = (1 - \xi) / 4 \\
 * \end{array}
 * @f]
 *
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 0 \qquad  \eta_{q0} = 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_quadrangle_4, _gt_quadrangle_4,
                                     _itp_lagrange_quadrangle_4, _ek_regular, 2,
                                     _git_segment, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2,
          aka::enable_if_t<aka::are_vectors<D1, D2>::value> *>
inline void InterpolationElement<_itp_lagrange_quadrangle_4>::computeShapes(
    const Eigen::MatrixBase<D1> &c, Eigen::MatrixBase<D2> &N) {
  N(0) = 1. / 4. * (1. - c(0)) * (1. - c(1)); /// N1(q_0)
  N(1) = 1. / 4. * (1. + c(0)) * (1. - c(1)); /// N2(q_0)
  N(2) = 1. / 4. * (1. + c(0)) * (1. + c(1)); /// N3(q_0)
  N(3) = 1. / 4. * (1. - c(0)) * (1. + c(1)); /// N4(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void InterpolationElement<_itp_lagrange_quadrangle_4>::computeDNDS(
    const Eigen::MatrixBase<D1> &c, Eigen::MatrixBase<D2> &dnds) {
  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial
   * \xi}
   *               & \frac{\partial N3}{\partial \xi}  & \frac{\partial
   * N4}{\partial \xi}\\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial
   * \eta}
   *               & \frac{\partial N3}{\partial \eta} & \frac{\partial
   * N4}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */

  dnds(0, 0) = -1. / 4. * (1. - c(1));
  dnds(0, 1) = 1. / 4. * (1. - c(1));
  dnds(0, 2) = 1. / 4. * (1. + c(1));
  dnds(0, 3) = -1. / 4. * (1. + c(1));

  dnds(1, 0) = -1. / 4. * (1. - c(0));
  dnds(1, 1) = -1. / 4. * (1. + c(0));
  dnds(1, 2) = 1. / 4. * (1. + c(0));
  dnds(1, 3) = 1. / 4. * (1. - c(0));
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_quadrangle_4>::computeD2NDS2(
    const vector_type & /*c*/, matrix_type &d2nds2) {
  d2nds2.zero();

  d2nds2(1, 0) = 1. / 4.;
  d2nds2(1, 1) = -1. / 4.;
  d2nds2(1, 2) = 1. / 4.;
  d2nds2(1, 3) = -1. / 4.;

  d2nds2(2, 0) = 1. / 4.;
  d2nds2(2, 1) = -1. / 4.;
  d2nds2(2, 2) = 1. / 4.;
  d2nds2(2, 3) = -1. / 4.;
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline Real GeometricalElement<_gt_quadrangle_4>::getInradius(
    const Eigen::MatrixBase<D> &coord) {
  auto &&u0 = coord.col(0);
  auto &&u1 = coord.col(1);
  auto &&u2 = coord.col(2);
  auto &&u3 = coord.col(3);
  Real a = (u0 - u1).norm();
  Real b = (u1 - u2).norm();
  Real c = (u2 - u3).norm();
  Real d = (u3 - u0).norm();

  // Real septimetre = (a + b + c + d) / 2.;

  // Real p = Math::distance_2d(coord + 0, coord + 4);
  // Real q = Math::distance_2d(coord + 2, coord + 6);

  // Real area = sqrt(4*(p*p * q*q) - (a*a + b*b + c*c + d*d)*(a*a + c*c - b*b -
  // d*d)) / 4.;
  // Real h = sqrt(area);  // to get a length
  // Real h = area / septimetre;  // formula of inradius for circumscritable
  // quadrelateral
  Real h = std::min({a, b, c, d});

  return h;
}
} // namespace akantu
