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
              q
   --x--------|--------x---> x
    -1        0        1
 @endverbatim
 *
 * @f{eqnarray*}{
 * w_1(x) &=& 1/2(1 - x) \\
 * w_2(x) &=& 1/2(1 + x)
 * @f}
 *
 * @f{eqnarray*}{
 * x_{q}  &=& 0
 * @f}
 */

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_segment_2, _gt_segment_2,
                                     _itp_lagrange_segment_2, _ek_regular, 1,
                                     _git_segment, 1);

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2,
          aka::enable_if_t<aka::are_vectors<D1, D2>::value> *>
inline void InterpolationElement<_itp_lagrange_segment_2>::computeShapes(
    const Eigen::MatrixBase<D1> &natural_coords, Eigen::MatrixBase<D2> &N) {

  /// natural coordinate
  Real c = natural_coords(0);
  /// shape functions
  N(0) = 0.5 * (1 - c);
  N(1) = 0.5 * (1 + c);
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void InterpolationElement<_itp_lagrange_segment_2>::computeDNDS(
    const Eigen::MatrixBase<D1> & /*natural_coords*/,
    Eigen::MatrixBase<D2> &dnds) {

  /// dN1/de
  dnds(0, 0) = -.5;
  /// dN2/de
  dnds(0, 1) = .5;
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_segment_2>::computeD2NDS2(
    const vector_type & /*natural_coords*/, matrix_type &d2nds2) {
  d2nds2.zero();
}

/* -------------------------------------------------------------------------- */
template <>
template <class D>
inline Real GeometricalElement<_gt_segment_2>::getInradius(
    const Eigen::MatrixBase<D> &coord) {
  return (coord.col(1) - coord.col(0)).norm();
}

/* -------------------------------------------------------------------------- */
template <>
template <class D1, class D2>
inline void
GeometricalElement<_gt_segment_2>::getNormal(const Eigen::MatrixBase<D1> &coord,
                                             Eigen::MatrixBase<D2> &normal) {
  assert(normal.size() == 2 && "The normal is only uniquely defined in 2D");
  Math::normal(coord.col(0) - coord.col(1), normal);
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
