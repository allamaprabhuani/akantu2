/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_GEOMETRY_UTILS_HH_
#define AKANTU_GEOMETRY_UTILS_HH_

namespace akantu {

class GeometryUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  static Vector<Real> outsideDirection(const Mesh & mesh,
                                       const Element & element);

  /// computes the normal on an element (assuming elements is flat)
  template <class Derived>
  static Vector<Real> normal(const Mesh & mesh,
                             const Eigen::MatrixBase<Derived> & coords,
                             const Element & element, bool outward = true);

  // computes normal at given covariant basis
  template <class Derived>
  static Vector<Real> normal(const Mesh & mesh, const Element & element,
                             Eigen::MatrixBase<Derived> & tangents,
                             bool outward = true);

  /// computes the orthogonal projection on a set of elements and
  /// returns natural projection and normal gap and index of element
  template <class Derived1, class Derived2, class Derived3, class Derived4,
            class ElementList>
  static Element orthogonalProjection(
      const Mesh & mesh, const Array<Real> & positions,
      const Eigen::MatrixBase<Derived1> & slave, const ElementList & elements,
      Real & gap, Eigen::MatrixBase<Derived2> & natural_projection,
      Eigen::MatrixBase<Derived3> & normal,
      Eigen::MatrixBase<Derived4> & tangent, Real alpha,
      Int max_iterations = 100, Real projection_tolerance = 1e-10,
      Real extension_tolerance = 1e-5);

  /// computes the natural projection on an element
  template <class Derived1, class Derived2>
  static std::pair<Vector<Real>, Vector<Real>>
  naturalProjection(const Eigen::MatrixBase<Derived1> & coords,
                    const Element & element,
                    const Eigen::MatrixBase<Derived2> & slave_coords,
                    Int max_iterations = 100, Real tolerance = 1e-10);

  /// computes the real projection on an element
  template <class Derived1, class Derived2, class Derived3>
  static Vector<Real>
  realProjection(const Eigen::MatrixBase<Derived1> & coords,
                 const Eigen::MatrixBase<Derived2> & slave,
                 const Eigen::MatrixBase<Derived3> & normal);

  /// computes the real projection from a natural coordinate
  template <class Derived1, class Derived2>
  static Vector<Real>
  realProjection(const Eigen::MatrixBase<Derived1> & coords,
                 const Element & element,
                 const Eigen::MatrixBase<Derived2> & natural_coord);

  /// computes the covariant basis/ local surface basis/ tangents on projection
  /// point
  template <class Derived1, class Derived2>
  static Matrix<Real>
  covariantBasis(const Eigen::MatrixBase<Derived1> & coords,
                 const Element & element,
                 Eigen::MatrixBase<Derived2> & natural_coord);

  /// computes the covariant basis/ local surface basis/ tangents on projection
  /// point
  template <class Derived1, class Derived2, class Derived3>
  static Matrix<Real>
  covariantBasis(const Eigen::MatrixBase<Derived1> & coords,
                 const Element & element,
                 const Eigen::MatrixBase<Derived2> & normal,
                 Eigen::MatrixBase<Derived3> & natural_coord);

  // computes the curvature on projection
  template <class Derived1, class Derived2>
  static Matrix<Real>
  curvature(const Eigen::MatrixBase<Derived1> & coords, const Element & element,
            const Eigen::MatrixBase<Derived2> & natural_coord);

  /// computes the contravariant basis on projection point
  template <class Derived>
  static Matrix<Real>
  contravariantBasis(const Eigen::MatrixBase<Derived> & covariant);

  /// computes metric tesnor with covariant components
  template <class Derived>
  static Matrix<Real>
  covariantMetricTensor(const Eigen::MatrixBase<Derived> & covariant_bases);

  /// computes metric tensor with contravariant components
  template <class Derived>
  static Matrix<Real>
  contravariantMetricTensor(const Eigen::MatrixBase<Derived> & covariant_bases);

  // computes curvature tensor with convariant components
  template <class Derived1, class Derived2, class Derived3>
  static Matrix<Real>
  covariantCurvatureTensor(const Eigen::MatrixBase<Derived1> & coords,
                           const Element & element,
                           const Eigen::MatrixBase<Derived2> & natural_coord,
                           const Eigen::MatrixBase<Derived3> & normal);

  /// checks if the element is truly a boundary element or not
  inline static bool isBoundaryElement(const Mesh & mesh,
                                       const Element & element);

  /// checks if the natural projection is valid for not
  template <class Derived>
  inline static bool
  isValidProjection(const Eigen::MatrixBase<Derived> & projection,
                    Real extension_tolerance = 1e-5);
};

} // namespace akantu

#include "geometry_utils_inline_impl.hh"

#endif /* AKANTU_GEOMETRY_UTILS_HH_ */
