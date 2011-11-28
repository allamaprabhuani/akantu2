/**
 * @file   aka_math.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 11:51:56 2010
 *
 * @brief  mathematical operations
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_MATH_H__
#define __AKANTU_AKA_MATH_H__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */

template<typename T>
class Vector;

class Math {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Matrix algebra                                                           */
  /* ------------------------------------------------------------------------ */
  /// @f$ y = A*x @f$
  static void matrix_vector(UInt m, UInt n,
			    const Vector<Real> & A,
			    const Vector<Real> & x,
			    Vector<Real> & y, Real alpha = 1.);

  /// @f$ y = A*x @f$
  static __aka_inline__ void matrix_vector(UInt m, UInt n,
				   const Real * A,
				   const Real * x,
				   Real * y, Real alpha = 1.);

  /// @f$ y = A^t*x @f$
  static __aka_inline__ void matrixt_vector(UInt m, UInt n,
				    const Real * A,
				    const Real * x,
				    Real * y, Real alpha = 1.);

  /// @f$ C = A*B @f$
  static void matrix_matrix(UInt m, UInt n, UInt k,
			   const Vector<Real> & A,
			   const Vector<Real> & B,
			   Vector<Real> & C, Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  static void matrix_matrixt(UInt m, UInt n, UInt k,
			     const Vector<Real> & A,
			     const Vector<Real> & B,
			     Vector<Real> & C, Real alpha = 1.);

  /// @f$ C = A*B @f$
  static __aka_inline__ void matrix_matrix(UInt m, UInt n, UInt k,
				   const Real * A,
				   const Real * B,
				   Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B @f$
  static __aka_inline__ void matrixt_matrix(UInt m, UInt n, UInt k,
				    const Real * A,
				    const Real * B,
				    Real * C, Real alpha = 1.);

  /// @f$ C = A*B^t @f$
  static __aka_inline__ void matrix_matrixt(UInt m, UInt n, UInt k,
				    const Real * A,
				    const Real * B,
				    Real * C, Real alpha = 1.);

  /// @f$ C = A^t*B^t @f$
  static __aka_inline__ void matrixt_matrixt(UInt m, UInt n, UInt k,
				     const Real * A,
				     const Real * B,
				     Real * C, Real alpha = 1.);

  template <bool tr_A, bool tr_B>
  static __aka_inline__ void matMul(UInt m, UInt n, UInt k,
			    Real alpha, const Real * A, const Real * B,
			    Real beta, Real * C);

  template <bool tr_A>
  static __aka_inline__ void matVectMul(UInt m, UInt n,
				Real alpha, const Real * A, const Real * x,
				Real beta, Real * y);

  static void matrix33_eigenvalues(Real * A,
				    Real * Adiag);

  /// determinent of a 3x3 matrix
  static __aka_inline__ Real det3(const Real * mat);

  /// determinent of a 2x2 matrix
  static __aka_inline__ Real det2(const Real * mat);

  /// inverse a 3x3 matrix
  static __aka_inline__ void inv3(const Real * mat, Real * inv);

  /// inverse a 2x2 matrix
  static __aka_inline__ void inv2(const Real * mat, Real * inv);

  /// vector cross product
  static __aka_inline__ void vectorProduct3(const Real * v1, const Real * v2, Real * res);

  /// compute normal a normal to a vector
  static __aka_inline__ void normal2(const Real * v1, Real * res);

  /// compute normal a normal to a vector
  static __aka_inline__ void normal3(const Real * v1,const Real * v2, Real * res);

  /// normalize a vector
  static __aka_inline__ void normalize2(Real * v);

  /// normalize a vector
  static __aka_inline__ void normalize3(Real * v);

  /// return norm of a 2-vector
  static __aka_inline__ Real norm2(const Real * v);

  /// return norm of a 3-vector
  static __aka_inline__ Real norm3(const Real * v);

  /// return the dot product between 2 vectors in 2d
  static __aka_inline__ Real vectorDot2(const Real * v1, const Real * v2);

  /// return the dot product between 2 vectors in 3d
  static __aka_inline__ Real vectorDot3(const Real * v1, const Real * v2);

  /* ------------------------------------------------------------------------ */
  /* Geometry                                                                 */
  /* ------------------------------------------------------------------------ */
  /// distance in 2D between x and y
  static __aka_inline__ Real distance_2d(const Real * x, const Real * y);

  /// distance in 3D between x and y
  static __aka_inline__ Real distance_3d(const Real * x, const Real * y);

  /// radius of the in-circle of a triangle
  static __aka_inline__ Real triangle_inradius(const Real * coord1, const Real * coord2, const Real * coord3);

  /// radius of the in-circle of a tetrahedron
  static __aka_inline__ Real tetrahedron_inradius(const Real * coord1, const Real * coord2, const Real * coord3, const Real * coord4);

  /// volume of a tetrahedron
  static __aka_inline__ Real tetrahedron_volume(const Real * coord1, const Real * coord2, const Real * coord3, const Real * coord4);

  /// compute the barycenter of n points
  static __aka_inline__ void barycenter(const Real * coord,
				UInt nb_points, UInt spatial_dimension,
				Real * barycenter);

  /// vector between x and y
  static __aka_inline__ void vector_2d(const Real * x, const Real * y, Real * vec);

  /// vector pointing from x to y in 3 spatial dimension
  static __aka_inline__ void vector_3d(const Real * x, const Real * y, Real * vec);

  /// test if two scalar are equal within a given tolerance
  static __aka_inline__ bool are_float_equal(Real x, Real y);

  /// test if two vectors are equal within a given tolerance
  static __aka_inline__ bool are_vector_equal(UInt n, Real * x, Real * y);

  /// test if a real is a NaN
  static __aka_inline__ bool isnan(Real x);

  /// test if the line x and y intersects each other
  static __aka_inline__ bool intersects(Real x_min, Real x_max, Real y_min, Real y_max);

  /// test if a is in the range [x_min, x_max]
  static __aka_inline__ bool is_in_range(Real a, Real x_min, Real x_max);

  static __aka_inline__ Real getTolerance() { return tolerance; };
  static __aka_inline__ void setTolerance(Real tol) { tolerance = tol; };

private:
  /// tolerance for functions that need one
  static Real tolerance;
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "aka_math_inline_impl.cc"
#endif


__END_AKANTU__

#endif /* __AKANTU_AKA_MATH_H__ */
