/**
 * @file   element_class_quadrangle_8_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed May 18 2011
 * @date last modification: Tue Sep 29 2020
 *
 * @brief  Specialization of the ElementClass for the _quadrangle_8
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
              \eta
               ^
               |
   (-1,1)    (0,1)   (1,1)
       x-------x-------x
       |       |       |
       |       |       |
       |       |       |
 (-1,0)|       |       |(1,0)
   ----x---------------X----->  \xi
       |       |       |
       |       |       |
       |       |       |
       |       |       |
       x-------x-------x
   (-1,-1)   (0,-1)  (1,-1)
               |
 @endverbatim
 *
 * @f[
 * \begin{array}{lll}
 * N1 = (1 - \xi) (1 - \eta)(- 1 - \xi - \eta) / 4
 *       & \frac{\partial N1}{\partial \xi}  = (1 - \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N1}{\partial \eta} = (1 - \xi)(\xi + 2 \eta) / 4 \\
 * N2 = (1 + \xi) (1 - \eta)(- 1 + \xi - \eta) / 4 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N2}{\partial \eta} = - (1 + \xi)(\xi - 2 \eta) / 4 \\
 * N3 = (1 + \xi) (1 + \eta)(- 1 + \xi + \eta) / 4 \\
 *       & \frac{\partial N3}{\partial \xi}  = (1 + \eta)(2 \xi + \eta) / 4
 *       & \frac{\partial N3}{\partial \eta} = (1 + \xi)(\xi + 2 \eta) / 4 \\
 * N4 = (1 - \xi) (1 + \eta)(- 1 - \xi + \eta) / 4
 *       & \frac{\partial N4}{\partial \xi}  = (1 + \eta)(2 \xi - \eta) / 4
 *       & \frac{\partial N4}{\partial \eta} = - (1 - \xi)(\xi - 2 \eta) / 4 \\
 * N5 = (1 - \xi^2) (1 - \eta) / 2
 *       & \frac{\partial N1}{\partial \xi}  = - \xi (1 - \eta)
 *       & \frac{\partial N1}{\partial \eta} = - (1 - \xi^2) / 2  \\
 * N6 = (1 + \xi) (1 - \eta^2) / 2 \\
 *       & \frac{\partial N2}{\partial \xi}  = (1 - \eta^2) / 2
 *       & \frac{\partial N2}{\partial \eta} = - \eta (1 + \xi) \\
 * N7 = (1 - \xi^2) (1 + \eta) / 2 \\
 *       & \frac{\partial N3}{\partial \xi}  = - \xi (1 + \eta)
 *       & \frac{\partial N3}{\partial \eta} = (1 - \xi^2) / 2 \\
 * N8 = (1 - \xi) (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \xi}  = - (1 - \eta^2) / 2
 *       & \frac{\partial N4}{\partial \eta} = - \eta (1 - \xi) \\
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
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_quadrangle_8, _gt_quadrangle_8,
                                     _itp_serendip_quadrangle_8, _ek_regular, 2,
                                     _git_segment, 3);

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2,
          aka::enable_if_t<aka::are_vectors<D1, D2>::value> *>
inline void InterpolationElement<_itp_serendip_quadrangle_8>::computeShapes(
    const Eigen::MatrixBase<D1> &c, Eigen::MatrixBase<D2> &N) {

  /// Natural coordinates
  const Real xi = c(0);
  const Real eta = c(1);

  N(0) = .25 * (1 - xi) * (1 - eta) * (-1 - xi - eta);
  N(1) = .25 * (1 + xi) * (1 - eta) * (-1 + xi - eta);
  N(2) = .25 * (1 + xi) * (1 + eta) * (-1 + xi + eta);
  N(3) = .25 * (1 - xi) * (1 + eta) * (-1 - xi + eta);
  N(4) = .5 * (1 - xi * xi) * (1 - eta);
  N(5) = .5 * (1 + xi) * (1 - eta * eta);
  N(6) = .5 * (1 - xi * xi) * (1 + eta);
  N(7) = .5 * (1 - xi) * (1 - eta * eta);
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void InterpolationElement<_itp_serendip_quadrangle_8>::computeDNDS(
    const Eigen::MatrixBase<D1> &c, Eigen::MatrixBase<D2> &dnds) {

  const Real xi = c(0);
  const Real eta = c(1);

  /// dN/dxi
  dnds(0, 0) = .25 * (1 - eta) * (2 * xi + eta);
  dnds(0, 1) = .25 * (1 - eta) * (2 * xi - eta);
  dnds(0, 2) = .25 * (1 + eta) * (2 * xi + eta);
  dnds(0, 3) = .25 * (1 + eta) * (2 * xi - eta);
  dnds(0, 4) = -xi * (1 - eta);
  dnds(0, 5) = .5 * (1 - eta * eta);
  dnds(0, 6) = -xi * (1 + eta);
  dnds(0, 7) = -.5 * (1 - eta * eta);

  /// dN/deta
  dnds(1, 0) = .25 * (1 - xi) * (2 * eta + xi);
  dnds(1, 1) = .25 * (1 + xi) * (2 * eta - xi);
  dnds(1, 2) = .25 * (1 + xi) * (2 * eta + xi);
  dnds(1, 3) = .25 * (1 - xi) * (2 * eta - xi);
  dnds(1, 4) = -.5 * (1 - xi * xi);
  dnds(1, 5) = -eta * (1 + xi);
  dnds(1, 6) = .5 * (1 - xi * xi);
  dnds(1, 7) = -eta * (1 - xi);
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
constexpr inline Real GeometricalElement<_gt_quadrangle_8>::getInradius(
    const Eigen::MatrixBase<D> &coord) {
  auto &&u0 = coord.col(0);
  auto &&u1 = coord.col(1);
  auto &&u2 = coord.col(2);
  auto &&u3 = coord.col(3);
  auto &&u4 = coord.col(4);
  auto &&u5 = coord.col(5);
  auto &&u6 = coord.col(6);
  auto &&u7 = coord.col(7);

  Real a = (u0 - u4).norm();
  Real b = (u4 - u1).norm();
  Real h = std::min(a, b);

  a = (u1 - u5).norm();
  b = (u5 - u2).norm();
  h = std::min(h, std::min(a, b));

  a = (u2 - u6).norm();
  b = (u6 - u3).norm();
  h = std::min(h, std::min(a, b));

  a = (u3 - u7).norm();
  b = (u7 - u0).norm();
  h = std::min(h, std::min(a, b));

  return h;
}

} // namespace akantu
