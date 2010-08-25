/**
 * @file   aka_math.h
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 11:51:56 2010
 *
 * @brief  mathematical operations
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_MATH_H__
#define __AKANTU_AKA_MATH_H__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Math {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* Matrix Algebra */
  /// y = A*x
  static void matrix_vector(UInt m, UInt n,
			    const Vector<Real> & A,
			    const Vector<Real> & x,
			    Vector<Real> & y);

  static inline void matrix_vector(UInt m, UInt n,
				   const Real * A,
				   const Real * x,
				   Real * y);

  /// C = A*B
  static void matrix_matrix(UInt m, UInt n, UInt k,
			   const Vector<Real> & A,
			   const Vector<Real> & B,
			   Vector<Real> & C);

  static inline void matrix_matrix(UInt m, UInt n, UInt k,
				   const Real * A,
				   const Real * B,
				   Real * C);

  /// C = A^t*B
  static inline void matrixt_matrix(UInt m, UInt n, UInt k,
				    const Real * A,
				    const Real * B,
				    Real * C);

  /// C = A*B^t
  static inline void matrix_matrixt(UInt m, UInt n, UInt k,
				    const Real * A,
				    const Real * B,
				    Real * C);

  /// C = A^t*B^t
  static inline void matrixt_matrixt(UInt m, UInt n, UInt k,
				     const Real * A,
				     const Real * B,
				     Real * C);

  /// determinent of a 3x3 matrix
  static inline Real det3(const Real * mat);

  /// determinent of a 2x2 matrix
  static inline Real det2(const Real * mat);

  /// inverse a 3x3 matrix
  static inline void inv3(const Real * mat, Real * inv);

  /// inverse a 2x2 matrix
  static inline void inv2(const Real * mat, Real * inv);

  /// vector cross product 
  static inline void vectorProduct3(const Real * v1, const Real * v2, Real * res);

  /// compute normal a normal to a vector 
  static inline void normal2(const Real * v1, Real * res);

  /// normalize a vector 
  static inline void normalize2(Real * v);

  /// normalize a vector 
  static inline Real norm2(Real * v);

  /* Geometry */
  /// distance in 2D between x and y
  static inline Real distance_2d(const Real * x, const Real * y);

  /// distance in 3D between x and y
  static inline Real distance_3d(const Real * x, const Real * y);

  /// radius of the in-circle of a triangle
  static inline Real triangle_inradius(const Real * coord);

  /// radius of the in-circle of a tetrahedron
  static inline Real tetrahedron_inradius(const Real * coord);

  /// volume of a tetrahedron
  static inline Real tetrahedron_volume(const Real * coord);

  ///
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "aka_math_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_AKA_MATH_H__ */
