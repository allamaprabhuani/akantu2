/**
 * @file   geometry_utils.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Sep 30 2019
 * @date last modification: Mon Sep 30 2019
 *
 * @brief  class to compute geometry related quantities 
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fe_engine.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GEOMETRY_UTILS_HH__
#define __AKANTU_GEOMETRY_UTILS_HH__

namespace akantu {

class GeometryUtils {
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// computes the normal on an element
  static void normal(const Mesh & mesh, const Array<Real> & positions,
		     const Element & element, Vector<Real> & normal, bool outward=true);

  /// computes the orthogonal projection on a set of elements and
  /// returns natural projection and normal gap and index of element
  static UInt orthogonalProjection(const Mesh & mesh, const Array<Real> & positions,
				   const Vector<Real> & slave,
				   const Array<Element> & elements,
				   Real & gap, Vector<Real> & natural_projection,
				   Vector<Real> & normal, Real alpha, Real tolerance = 1e-10);

  /// computes the natural projection on an element
  static void naturalProjection(const Mesh & mesh, const Array<Real> & positions,
				const Element & element, Vector<Real> & real_projection,
				Vector<Real> & natural_projection, Real tolerance = 1e-10);
  
  /// computes the real projection on an element
  static void realProjection(const Mesh & mesh, const Array<Real> & positions,
				    const Vector<Real> & slave,  const Element & element,
				    const Vector<Real> & normal, Vector<Real> & projection);

  /// computes the real projection from a natural coordinate
  static void realProjection(const Mesh & mesh, const Array<Real> & positions,
			     const Element & element, const Vector<Real> & natural_coord,
			     Vector<Real> & projection);
  
  /// computes the covariant basis/ local surface basis/ tangents on projection
  /// point
  static void covariantBasis(const Mesh & mesh, const Array<Real> & positions,
			     const Element & element,  const Vector<Real> & normal,
			     Vector<Real> & natural_coord,
			     Matrix<Real> & basis);

  // computes the curvature on projection
  static void curvature(const Mesh & mesh, const Array<Real> & positions,
			const Element & element, const Vector<Real> & natural_coord,
			Matrix<Real> & curvature);

  
  /// computes the contravariant basis on projection point
  static void contravariantBasis(const Matrix<Real> & covariant,
				 Matrix<Real> & contravariant);

  /// computes metric tesnor with covariant components
  static Matrix<Real> covariantMetricTensor(const Matrix<Real> & );

  /// computes metric tensor with contravariant components
  static Matrix<Real> contravariantMetricTensor(const Matrix<Real> & );

  // computes curvature tensor with convariant components
  static Matrix<Real> covariantCurvatureTensor(const Mesh &,
					       const Array<Real> &,
					       const Element &,
					       const Vector<Real> &,
					       const Vector<Real> &);  
  
  /// checks if the element is truly a boundary element or not
  inline static bool isBoundaryElement(const Mesh & mesh, const Element & element);

  /// checks if the natural projection is valid for not
  inline static bool isValidProjection(const Vector<Real> & projection);


};
  
}


#include "geometry_utils_inline_impl.cc"  


#endif /* __AKANTU_GEOMETRY_UTILS_HH__ */
