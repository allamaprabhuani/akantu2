/**
 * @file   element_class_hermite_inline_impl.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Nov 10 2017
 * @date last modification: Tue Feb 09 2021
 *
 * @brief  Specialization of the element_class class for the type
 * _hermite
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * @f[
 *   \begin{array}{ll}
 *     M_1(\xi) &= 1/4(\xi^{3}/-3\xi+2)\\
 *     M_2(\xi) &= -1/4(\xi^{3}-3\xi-2)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L_1(\xi) &= 1/4(\xi^{3}-\xi^{2}-\xi+1)\\
 *     L_2(\xi) &= 1/4(\xi^{3}+\xi^{2}-\xi-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     M'_1(\xi) &= 3/4(\xi^{2}-1)\\
 *     M'_2(\xi) &= -3/4(\xi^{2}-1)
 *   \end{array}
 *
 *   \begin{array}{ll}
 *     L'_1(\xi) &= 1/4(3\xi^{2}-2\xi-1)\\
 *     L'_2(\xi) &= 1/4(3\xi^{2}+2\xi-1)
 *   \end{array}
 *@f]
 *
 *
 *@f[
 * \begin{array}{ll}
 *   N'_1(\xi) &= -1/2\\
 *   N'_2(\xi) &= 1/2
 * \end{array}]
 *
 * \begin{array}{ll}
 *   -M''_1(\xi) &= -3\xi/2\\
 *   -M''_2(\xi) &= 3\xi/2\\
 * \end{array}
 *
 * \begin{array}{ll}
 *   -L''_1(\xi) &= -1/2a(3\xi/a-1)\\
 *   -L''_2(\xi) &= -1/2a(3\xi/a+1)
 * \end{array}
 *@f]
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_static_if.hh"
//#include "element_class_structural.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_HH_
#define AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */

AKANTU_DEFINE_STRUCTURAL_INTERPOLATION_TYPE_PROPERTY(_itp_hermite_2,
                                                     _itp_lagrange_segment_2, 2,
                                                     1, 4);

/* -------------------------------------------------------------------------- */
namespace details {
  template <class D1>
  inline Real computeLength(const Eigen::MatrixBase<D1> & real_coord) {
    auto && x1 = real_coord(0);
    auto && x2 = real_coord(1);
    return x1.distance(x2);
  }

  template <class D1, class D2>
  inline void computeShapes(const Eigen::MatrixBase<D1> & natural_coords, Real a,
                            Eigen::MatrixBase<D2> & N) {
    /// natural coordinate
    Real xi = natural_coords(0);
    auto xi2 = xi * xi;
    auto xi3 = xi * xi * xi;
    // Cubic Hermite splines interpolating displacement
    auto M1 = 1. / 4. * (2. - 3. * xi + xi3);
    auto M2 = 1. / 4. * (2. + 3. * xi - xi3);
    auto L1 = a / 4. * (1 - xi - xi2 + xi3);
    auto L2 = a / 4. * (-1 - xi + xi2 + xi3);

#if 1 // Version where we also interpolate the rotations
      // Derivatives (with respect to x) of previous functions interpolating
      // rotations
    auto M1_ = 3. / (4. * a) * (xi2 - 1);
    auto M2_ = 3. / (4. * a) * (1 - xi2);
    auto L1_ = 1 / 4. * (3 * xi2 - 2 * xi - 1);
    auto L2_ = 1 / 4. * (3 * xi2 + 2 * xi - 1);

    // clang-format off
      //    v1   t1   v2   t2
      N << M1 , L1 , M2 , L2,   // displacement interpolation
           M1_, L1_, M2_, L2_; // rotation interpolation
    // clang-format on

#else // Version where we only interpolate displacements
      // clang-format off
      //    v1  t1  v2  t2
      N = {{M1, L1, M2, L2}};
// clang-format on
#endif
  }

  /* ---------------------------------------------------------------------- */
  template <class D1, class D2>
  inline void computeDNDS(const Eigen::MatrixBase<D1> & natural_coords, Real a,
                          Eigen::MatrixBase<D2> & B) {
    // natural coordinate
    Real xi = natural_coords(0);
    // Derivatives with respect to xi for rotations
    auto M1 = 3. / 2. * xi;
    auto M2 = 3. / 2. * (-xi);
    auto L1 = 1. * a / 2. * (3 * xi - 1);
    auto L2 = 1. * a / 2. * (3 * xi + 1);

    //   v1  t1  v2  t2
    B << M1, L1, M2, L2; // computing curvature : {chi} = [B]{d}
    B /= a;              // to account for first order deriv w/r to x
  }
} // namespace details

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
inline void
InterpolationElement<_itp_hermite_2, _itk_structural>::computeShapes(
    const Eigen::MatrixBase<D1> & natural_coords,
    const Eigen::MatrixBase<D2> & real_coord, Eigen::MatrixBase<D3> & N) {
  auto L = details::computeLength(real_coord);
  details::computeShapes(natural_coords, L / 2, N);
}

/* -------------------------------------------------------------------------- */
template <>
template <typename D1, typename D2, typename D3>
inline void InterpolationElement<_itp_hermite_2, _itk_structural>::computeDNDS(
    const Eigen::MatrixBase<D1> & Xs, const Eigen::MatrixBase<D2> & xs,
    Eigen::MatrixBase<D3> & B) {
  auto L = details::computeLength(xs);
  details::computeDNDS(Xs, L / 2, B);
}

} // namespace akantu
#endif /* AKANTU_ELEMENT_CLASS_HERMITE_INLINE_IMPL_HH_ */
