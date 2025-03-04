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
      x
    (0)
 @endverbatim
 *
 * @f{eqnarray*}{
 * N1 &=& 1
 * @f}
 *
 * @f{eqnarray*}{
 * q_0 &=& 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_point_1, _gt_point, _itp_lagrange_point_1,
                                     _ek_regular, 0, _git_point, 1);

template <>
template <class D1, class D2, class D3>
inline void
ElementClass<_point_1, _ek_regular>::computeNormalsOnNaturalCoordinates(
    const Eigen::MatrixBase<D1> & /*coord*/,
    const Eigen::MatrixBase<D2> & /*f*/, Eigen::MatrixBase<D3> & /*normals*/) {}

/* --------------r------------------------------------------------------------
 */
template <>
template <class D1, class D2,
          aka::enable_if_t<aka::are_vectors<D1, D2>::value> *>
inline void InterpolationElement<_itp_lagrange_point_1>::computeShapes(
    const Eigen::MatrixBase<D1> & /*natural_coords*/,
    Eigen::MatrixBase<D2> & N) {
  N(0) = 1; /// N1(q_0)
}
/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void InterpolationElement<_itp_lagrange_point_1>::computeDNDS(
    const Eigen::MatrixBase<D1> & /*natural_coords*/,
    Eigen::MatrixBase<D2> & /*dnds*/) {}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline Real InterpolationElement<_itp_lagrange_point_1>::computeSpecialJacobian(
    const Eigen::MatrixBase<D> & /*J*/) {
  return 0.;
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline Real GeometricalElement<_gt_point>::getInradius(
    const Eigen::MatrixBase<D> & /*coord*/) {
  return 0.;
}
} // namespace akantu
