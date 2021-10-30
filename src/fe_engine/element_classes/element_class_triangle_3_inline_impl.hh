/**
 * @file   element_class_triangle_3_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Specialization of the element_class class for the type _triangle_3
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
     x (0,0,1)
     |`
     |  `
     |  q `
     |  °   `
     x--------x----->  \xi
    (1,0,0)      (0,1,0)
 @endverbatim
 *
 * @f{eqnarray*}{
 * N1 &=& 1 - \xi - \eta \\
 * N2 &=& \xi \\
 * N3 &=& \eta
 * @f}
 *
 * @f{eqnarray*}{
 * \xi_{q0}  &=& 1/3 \qquad  \eta_{q0} = 1/3
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_triangle_3, _gt_triangle_3,
                                     _itp_lagrange_triangle_3, _ek_regular, 2,
                                     _git_triangle, 1);

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2,
          aka::enable_if_t<aka::are_vectors<D1, D2>::value> *>
inline void InterpolationElement<_itp_lagrange_triangle_3>::computeShapes(
    const Eigen::MatrixBase<D1> & X, Eigen::MatrixBase<D2> & N) {
  N(0) = 1 - X(0) - X(1); /// @f$ c0 = 1 - \xi - \eta @f$
  N(1) = X(0);            /// @f$ c1 = \xi @f$
  N(2) = X(1);            /// @f$ c2 = \eta @f$
}
/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void InterpolationElement<_itp_lagrange_triangle_3>::computeDNDS(
    const Eigen::MatrixBase<D1> & /*natural_coords*/,
    Eigen::MatrixBase<D2> & dnds) {
  /**
   * @f[
   * dnds = \left(
   *          \begin{array}{cccccc}
   *            \frac{\partial N1}{\partial \xi}  & \frac{\partial N2}{\partial
   * \xi}  & \frac{\partial N3}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta} & \frac{\partial N2}{\partial
   * \eta} & \frac{\partial N3}{\partial \eta}
   *          \end{array}
   *        \right)
   * @f]
   */
  dnds(0, 0) = -1.;
  dnds(0, 1) = 1.;
  dnds(0, 2) = 0.;
  dnds(1, 0) = -1.;
  dnds(1, 1) = 0.;
  dnds(1, 2) = 1.;
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_triangle_3>::computeD2NDS2(
    const vector_type & /*natural_coords*/, matrix_type & d2nds2) {
  d2nds2.zero();
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline Real GeometricalElement<_gt_triangle_3>::getInradius(
    const Eigen::MatrixBase<D> & coord) {
  auto && coord1 = coord.col(0);
  auto && coord2 = coord.col(1);
  auto && coord3 = coord.col(2);

  return Math::triangle_inradius(coord1, coord2, coord3);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
