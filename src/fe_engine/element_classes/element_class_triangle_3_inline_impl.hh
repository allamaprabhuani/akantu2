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
inline void InterpolationElement<_itp_lagrange_triangle_3>::computeShapes(
    const Ref<const VectorXr> & natural_coords, Ref<VectorXr> N) {

  /// Natural coordinates
  Real c0 =
      1 - natural_coords(0) - natural_coords(1); /// @f$ c0 = 1 - \xi - \eta @f$
  Real c1 = natural_coords(0);                   /// @f$ c1 = \xi @f$
  Real c2 = natural_coords(1);                   /// @f$ c2 = \eta @f$

  N(0) = c0; /// N1(q_0)
  N(1) = c1; /// N2(q_0)
  N(2) = c2; /// N3(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
inline void InterpolationElement<_itp_lagrange_triangle_3>::computeDNDS(
    const Ref<const VectorXr> & /*natural_coords*/, Ref<MatrixXr> dnds) {
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
inline void
InterpolationElement<_itp_lagrange_triangle_3>::computeSpecialJacobian(
    const Ref<const MatrixXr> & J, Real & jac) {
  Eigen::Map<const Eigen::Matrix<Real, 2, 3>> Jstatic(
      J.data());
  jac = Jstatic.row(0).cross(Jstatic.row(1)).norm();
}

/* -------------------------------------------------------------------------- */
template <>
inline Real GeometricalElement<_gt_triangle_3>::getInradius(
    const Ref<const MatrixXr> & coord) {
  return 2. * Math::triangle_inradius(coord.col(0), coord.col(1),
                                      coord.col(2));
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
