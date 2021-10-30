/**
 * @file   aka_math_tmpl.hh
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Fri Dec 11 2020
 *
 * @brief  Implementation of the inline functions of the math toolkit
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
#include "aka_blas_lapack.hh"
#include "aka_math.hh"
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <typeinfo>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace Math {
  /* ------------------------------------------------------------------------ */
  template <typename T>
  inline void solve(UInt n, const T * A, T * x, const T * b) {
    int N = n;
    int info;
    auto * ipiv = new int[N];
    auto * lu_A = new T[N * N];

    std::copy(A, A + N * N, lu_A);

    aka_getrf(&N, &N, lu_A, &N, ipiv, &info);
    if (info > 0) {
      AKANTU_ERROR("Singular matrix - cannot factorize it (info: " << info
                                                                   << " )");
    }

    char trans = 'N';
    int nrhs = 1;

    std::copy(b, b + N, x);

    aka_getrs(&trans, &N, &nrhs, lu_A, &N, ipiv, x, &N, &info);
    if (info != 0) {
      AKANTU_ERROR("Cannot solve the system (info: " << info << " )");
    }

    delete[] ipiv;
    delete[] lu_A;
  }

  /* ------------------------------------------------------------------------ */
  inline void vectorProduct3(const Real * v1, const Real * v2, Real * res) {
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  /* ------------------------------------------------------------------------ */
  inline Real vectorDot2(const Real * v1, const Real * v2) {
    return (v1[0] * v2[0] + v1[1] * v2[1]);
  }

  /* ------------------------------------------------------------------------ */
  inline Real vectorDot3(const Real * v1, const Real * v2) {
    return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
  }

  /* ------------------------------------------------------------------------ */
  inline Real distance_2d(const Real * x, const Real * y) {
    return std::sqrt((y[0] - x[0]) * (y[0] - x[0]) +
                     (y[1] - x[1]) * (y[1] - x[1]));
  }

  /* ------------------------------------------------------------------------ */
  template <class D1, class D2, class D3>
  static inline Real triangle_inradius(const Eigen::MatrixBase<D1> & coord1,
                                       const Eigen::MatrixBase<D2> & coord2,
                                       const Eigen::MatrixBase<D3> & coord3) {
    auto a = coord1.distance(coord2);
    auto b = coord2.distance(coord3);
    auto c = coord1.distance(coord3);

    auto s = (a + b + c) * 0.5;

    return std::sqrt((s - a) * (s - b) * (s - c) / s);
  }

  /* ------------------------------------------------------------------------ */
  template <class D1, class D2, class D3, class D4>
  inline Real tetrahedron_volume(const Eigen::MatrixBase<D1> & coord1,
                                 const Eigen::MatrixBase<D2> & coord2,
                                 const Eigen::MatrixBase<D3> & coord3,
                                 const Eigen::MatrixBase<D4> & coord4) {
    Matrix<Real, 3, 3> xx;

    xx.col(0) = coord2;
    xx.col(1) = coord3;
    xx.col(2) = coord4;
    auto vol = xx.determinant();

    xx.col(0) = coord1;
    vol -= xx.determinant();

    xx.col(1) = coord2;
    vol += xx.determinant();

    xx.col(2) = coord3;
    vol -= xx.determinant();

    vol /= 6;

    return vol;
  }

  /* ------------------------------------------------------------------------ */
  template <class D1, class D2, class D3, class D4>
  inline Real tetrahedron_inradius(const Eigen::MatrixBase<D1> & coord1,
                                   const Eigen::MatrixBase<D2> & coord2,
                                   const Eigen::MatrixBase<D3> & coord3,
                                   const Eigen::MatrixBase<D4> & coord4) {
    auto l12 = coord1.distance(coord2);
    auto l13 = coord1.distance(coord3);
    auto l14 = coord1.distance(coord4);
    auto l23 = coord2.distance(coord3);
    auto l24 = coord2.distance(coord4);
    auto l34 = coord3.distance(coord4);

    auto s1 = (l12 + l23 + l13) * 0.5;
    s1 = std::sqrt(s1 * (s1 - l12) * (s1 - l23) * (s1 - l13));

    auto s2 = (l12 + l24 + l14) * 0.5;
    s2 = std::sqrt(s2 * (s2 - l12) * (s2 - l24) * (s2 - l14));

    auto s3 = (l23 + l34 + l24) * 0.5;
    s3 = std::sqrt(s3 * (s3 - l23) * (s3 - l34) * (s3 - l24));

    auto s4 = (l13 + l34 + l14) * 0.5;
    s4 = std::sqrt(s4 * (s4 - l13) * (s4 - l34) * (s4 - l14));

    auto volume = Math::tetrahedron_volume(coord1, coord2, coord3, coord4);

    return 3 * volume / (s1 + s2 + s3 + s4);
  }

  /* ------------------------------------------------------------------------ */
  template <class D1, class D2>
  inline void barycenter(const Eigen::MatrixBase<D1> &coord,
                         Eigen::MatrixBase<D2> &barycenter) {
    barycenter.zero();
    for (auto && x : coord) {
      barycenter += x;
    }
    barycenter /= (Real)coord.cols();
  }

  /* ------------------------------------------------------------------------ */
  inline void vector_2d(const Real * x, const Real * y, Real * res) {
    res[0] = y[0] - x[0];
    res[1] = y[1] - x[1];
  }

  /* ------------------------------------------------------------------------ */
  inline void vector_3d(const Real * x, const Real * y, Real * res) {
    res[0] = y[0] - x[0];
    res[1] = y[1] - x[1];
    res[2] = y[2] - x[2];
  }

  /* ------------------------------------------------------------------------ */
  /// Combined absolute and relative tolerance test proposed in
  /// Real-time collision detection by C. Ericson (2004)
  inline bool are_float_equal(const Real x, const Real y) {
    Real abs_max = std::max(std::abs(x), std::abs(y));
    abs_max = std::max(abs_max, Real(1.));
    return std::abs(x - y) <= (tolerance * abs_max);
  }

  /* ------------------------------------------------------------------------ */
  inline bool isnan(Real x) {
#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable : 1572)
#endif // defined(__INTEL_COMPILER)

    // x = x return false means x = quiet_NaN
    return !(x == x);

#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif // defined(__INTEL_COMPILER)
  }

  /* ------------------------------------------------------------------------ */
  inline bool are_vector_equal(UInt n, Real * x, Real * y) {
    bool test = true;
    for (Int i = 0; i < n; ++i) {
      test &= are_float_equal(x[i], y[i]);
    }

    return test;
  }

  /* ------------------------------------------------------------------------ */
  inline bool intersects(Real x_min, Real x_max, Real y_min, Real y_max) {
    return not((x_max < y_min) or (x_min > y_max));
  }

  /* ------------------------------------------------------------------------ */
  inline bool is_in_range(Real a, Real x_min, Real x_max) {
    return ((a >= x_min) and (a <= x_max));
  }

  /* ------------------------------------------------------------------------ */
  template <UInt p, typename T> inline T pow(T x) {
    return (pow<p - 1, T>(x) * x);
  }
  template <> inline UInt pow<0, UInt>(__attribute__((unused)) UInt x) {
    return (1);
  }
  template <> inline Real pow<0, Real>(__attribute__((unused)) Real x) {
    return (1.);
  }

  /* ------------------------------------------------------------------------ */

  template <class Functor>
  Real NewtonRaphson::solve(const Functor & funct, Real x_0) {
    Real x = x_0;
    Real f_x = funct.f(x);
    UInt iter = 0;
    while (std::abs(f_x) > this->tolerance && iter < this->max_iteration) {
      x -= f_x / funct.f_prime(x);
      f_x = funct.f(x);
      iter++;
    }

    AKANTU_DEBUG_ASSERT(iter < this->max_iteration,
                        "Newton Raphson ("
                            << funct.name << ") solve did not converge in "
                            << this->max_iteration << " iterations (tolerance: "
                            << this->tolerance << ")");

    return x;
  }

} // namespace Math
} // namespace akantu
