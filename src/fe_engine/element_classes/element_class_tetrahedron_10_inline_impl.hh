/**
 * @file   element_class_tetrahedron_10_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 *
 * @date creation: Fri Jul 16 2010
 * @date last modification: Fri Feb 07 2020
 *
 * @brief  Specialization of the element_class class for the type
 * _tetrahedron_10
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

/* -------------------------------------------------------------------------- */
#include "element_class.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
AKANTU_DEFINE_ELEMENT_CLASS_PROPERTY(_tetrahedron_10, _gt_tetrahedron_10,
                                     _itp_lagrange_tetrahedron_10, _ek_regular,
                                     3, _git_tetrahedron, 2);

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline void InterpolationElement<_itp_lagrange_tetrahedron_10>::computeShapes(
    const vector_type & natural_coords, vector_type & N) {
  /// Natural coordinates
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);
  Real zeta = natural_coords(2);
  Real sum = xi + eta + zeta;
  Real c0 = 1 - sum;
  Real c1 = 1 - 2 * sum;
  Real c2 = 2 * xi - 1;
  Real c3 = 2 * eta - 1;
  Real c4 = 2 * zeta - 1;

  /// Shape functions
  N(0) = c0 * c1;
  N(1) = xi * c2;
  N(2) = eta * c3;
  N(3) = zeta * c4;
  N(4) = 4 * xi * c0;
  N(5) = 4 * xi * eta;
  N(6) = 4 * eta * c0;
  N(7) = 4 * zeta * c0;
  N(8) = 4 * xi * zeta;
  N(9) = 4 * eta * zeta;
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type, class matrix_type>
inline void InterpolationElement<_itp_lagrange_tetrahedron_10>::computeDNDS(
    const vector_type & natural_coords, matrix_type & dnds) {
  /**
   * \f[
   * dnds = \left(
   *          \begin{array}{cccccccccc}
   *            \frac{\partial N1}{\partial \xi}   & \frac{\partial N2}{\partial
   * \xi}
   *          & \frac{\partial N3}{\partial \xi}   & \frac{\partial N4}{\partial
   * \xi}
   *          & \frac{\partial N5}{\partial \xi}   & \frac{\partial N6}{\partial
   * \xi}
   *          & \frac{\partial N7}{\partial \xi}   & \frac{\partial N8}{\partial
   * \xi}
   *          & \frac{\partial N9}{\partial \xi}   & \frac{\partial
   * N10}{\partial \xi} \\
   *            \frac{\partial N1}{\partial \eta}  & \frac{\partial N2}{\partial
   * \eta}
   *          & \frac{\partial N3}{\partial \eta}  & \frac{\partial N4}{\partial
   * \eta}
   *          & \frac{\partial N5}{\partial \eta}  & \frac{\partial N6}{\partial
   * \eta}
   *          & \frac{\partial N7}{\partial \eta}  & \frac{\partial N8}{\partial
   * \eta}
   *          & \frac{\partial N9}{\partial \eta}  & \frac{\partial
   * N10}{\partial \eta} \\
   *            \frac{\partial N1}{\partial \zeta} & \frac{\partial N2}{\partial
   * \zeta}
   *          & \frac{\partial N3}{\partial \zeta} & \frac{\partial N4}{\partial
   * \zeta}
   *          & \frac{\partial N5}{\partial \zeta} & \frac{\partial N6}{\partial
   * \zeta}
   *          & \frac{\partial N7}{\partial \zeta} & \frac{\partial N8}{\partial
   * \zeta}
   *          & \frac{\partial N9}{\partial \zeta} & \frac{\partial
   * N10}{\partial \zeta}
   *          \end{array}
   *        \right)
   * \f]
   */

  /// Natural coordinates
  Real xi = natural_coords(0);
  Real eta = natural_coords(1);
  Real zeta = natural_coords(2);
  Real sum = xi + eta + zeta;

  /// \frac{\partial N_i}{\partial \xi}
  dnds(0, 0) = 4 * sum - 3;
  dnds(0, 1) = 4 * xi - 1;
  dnds(0, 2) = 0;
  dnds(0, 3) = 0;
  dnds(0, 4) = 4 * (1 - sum - xi);
  dnds(0, 5) = 4 * eta;
  dnds(0, 6) = -4 * eta;
  dnds(0, 7) = -4 * zeta;
  dnds(0, 8) = 4 * zeta;
  dnds(0, 9) = 0;

  /// \frac{\partial N_i}{\partial \eta}
  dnds(1, 0) = 4 * sum - 3;
  dnds(1, 1) = 0;
  dnds(1, 2) = 4 * eta - 1;
  dnds(1, 3) = 0;
  dnds(1, 4) = -4 * xi;
  dnds(1, 5) = 4 * xi;
  dnds(1, 6) = 4 * (1 - sum - eta);
  dnds(1, 7) = -4 * zeta;
  dnds(1, 8) = 0;
  dnds(1, 9) = 4 * zeta;

  /// \frac{\partial N_i}{\partial \zeta}
  dnds(2, 0) = 4 * sum - 3;
  dnds(2, 1) = 0;
  dnds(2, 2) = 0;
  dnds(2, 3) = 4 * zeta - 1;
  dnds(2, 4) = -4 * xi;
  dnds(2, 5) = 0;
  dnds(2, 6) = -4 * eta;
  dnds(2, 7) = 4 * (1 - sum - zeta);
  dnds(2, 8) = 4 * xi;
  dnds(2, 9) = 4 * eta;
}

/* -------------------------------------------------------------------------- */
template <>
inline Real GeometricalElement<_gt_tetrahedron_10>::getInradius(
    const Matrix<Real> & coord) {
  // Only take the four corner tetrahedra
  UInt tetrahedra[4][4] = {
      {0, 4, 6, 7}, {4, 1, 5, 8}, {6, 5, 2, 9}, {7, 8, 9, 3}};

  Real inradius = std::numeric_limits<Real>::max();
  for (UInt t = 0; t < 4; t++) {
    Real ir = Math::tetrahedron_inradius(
        coord(tetrahedra[t][0]).storage(), coord(tetrahedra[t][1]).storage(),
        coord(tetrahedra[t][2]).storage(), coord(tetrahedra[t][3]).storage());
    inradius = std::min(ir, inradius);
  }

  return 2. * inradius;
}
} // namespace akantu
