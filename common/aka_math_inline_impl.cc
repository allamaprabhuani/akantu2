/**
 * @file   aka_math_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 13:20:35 2010
 *
 * @brief  Implementation of the inline functions of the math toolkit
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */
#ifdef AKANTU_USE_CBLAS
# ifndef AKANTU_USE_CBLAS_MKL
#  include <cblas.h>
# else // AKANTU_USE_CBLAS_MKL
#  include <mkl_cblas.h>
# endif //AKANTU_USE_CBLAS_MKL
#endif //AKANTU_USE_CBLAS

/* -------------------------------------------------------------------------- */
inline void Math::matrix_vector(UInt m, UInt n,
				const Real * A,
				const Real * x,
				Real * y) {
#ifdef AKANTU_USE_CBLAS
  /// y = alpha*op(A)*x + beta*y
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      m, n, 1, A, n, x, 1, 0, y, 1);
#else
  memset(y, 0, m*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    UInt A_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      y[i] += A[A_i + j] * x[j];
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrix(UInt m, UInt n, UInt k,
				const Real * A,
				const Real * B,
				Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      m, n, k,
	      1,
	      A, k,
	      B, n,
	      0,
	      C, n);
#else
  memset(C, 0, m*n*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    UInt A_i = i * k;
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[A_i + l] * B[l * n + j];
      }
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrix(UInt m, UInt n, UInt k,
				 const Real * A,
				 const Real * B,
				 Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	      m, n, k,
	      1,
	      A, m,
	      B, n,
	      0,
	      C, n);
#else
  memset(C, 0, m*n*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    //    UInt A_i = i * k;
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[l * m + i] * B[l * n + j];
      }
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrixt(UInt m, UInt n, UInt k,
				const Real * A,
				const Real * B,
				Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	      m, n, k,
	      1,
	      A, k,
	      B, k,
	      0,
	      C, n);
#else
  memset(C, 0, m*n*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    UInt A_i = i * k;
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      UInt B_j = j * k;
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[A_i + l] * B[B_j + l];
      }
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrixt(UInt m, UInt n, UInt k,
				  const Real * A,
				  const Real * B,
				  Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	      m, n, k,
	      1,
	      A, m,
	      B, k,
	      0,
	      C, n);
#else
  memset(C, 0, m * n * sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      UInt B_j = j * k;
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[l * m + i] * B[B_j + l];
      }
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline Real Math::det2(const Real * mat) {
  return mat[0]*mat[3] - mat[1]*mat[2];
}

/* -------------------------------------------------------------------------- */
inline Real Math::det3(const Real * mat) {
  return
      mat[0]*(mat[4]*mat[8]-mat[7]*mat[5])
    - mat[3]*(mat[1]*mat[8]-mat[7]*mat[2])
    + mat[6]*(mat[1]*mat[5]-mat[4]*mat[2]);
}

/* -------------------------------------------------------------------------- */
inline void Math::normal2(const Real * vec,Real * normal) {
    normal[0] = vec[1];
    normal[1] = -vec[0];
    Math::normalize2(normal);
}

/* -------------------------------------------------------------------------- */
inline void Math::normal3(const Real * vec1,const Real * vec2,Real * normal) {
  Math::vectorProduct3(vec1,vec2,normal);
  Math::normalize3(normal);
}

/* -------------------------------------------------------------------------- */
inline void Math::normalize2(Real * vec) {
  Real norm = Math::norm2(vec);
  vec[0] /= norm;
  vec[1] /= norm;
}

/* -------------------------------------------------------------------------- */
inline void Math::normalize3(Real * vec) {
  Real norm = Math::norm3(vec);
  vec[0] /= norm;
  vec[1] /= norm;
  vec[2] /= norm;
}

/* -------------------------------------------------------------------------- */
inline Real Math::norm2(const Real * vec) {
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::norm3(const Real * vec) {
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

/* -------------------------------------------------------------------------- */
inline void Math::inv2(const Real * mat,Real * inv) {
  Real det_mat = det2(mat);

  inv[0] =  mat[3] / det_mat;
  inv[1] = -mat[1] / det_mat;
  inv[2] = -mat[2] / det_mat;
  inv[3] =  mat[0] / det_mat;
}

/* -------------------------------------------------------------------------- */
inline void Math::inv3(const Real * mat,Real * inv) {
  Real det_mat = det3(mat);

  inv[0] = (mat[4]*mat[8] - mat[7]*mat[5])/det_mat;
  inv[1] = (mat[2]*mat[7] - mat[8]*mat[1])/det_mat;
  inv[2] = (mat[1]*mat[5] - mat[4]*mat[2])/det_mat;
  inv[3] = (mat[5]*mat[6] - mat[8]*mat[3])/det_mat;
  inv[4] = (mat[0]*mat[8] - mat[6]*mat[2])/det_mat;
  inv[5] = (mat[2]*mat[3] - mat[5]*mat[0])/det_mat;
  inv[6] = (mat[3]*mat[7] - mat[6]*mat[4])/det_mat;
  inv[7] = (mat[1]*mat[6] - mat[7]*mat[0])/det_mat;
  inv[8] = (mat[0]*mat[4] - mat[3]*mat[1])/det_mat;
}

/* -------------------------------------------------------------------------- */
inline void Math::vectorProduct3(const Real * v1, const Real * v2, Real * res) {
  res[0] = v1[1]*v2[2] - v1[2]*v2[1];
  res[1] = v1[2]*v2[0] - v1[0]*v2[2];
  res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* -------------------------------------------------------------------------- */
inline Real Math::distance_2d(const Real * x, const Real * y) {
  return sqrt((y[0] - x[0])*(y[0] - x[0]) + (y[1] - x[1])*(y[1] - x[1]));
}

/* -------------------------------------------------------------------------- */
inline Real Math::triangle_inradius(const Real * coord1,
				    const Real * coord2,
				    const Real * coord3) {
  /**
   * @f{eqnarray*}{
   * r &=& A / s \\
   * A &=& 1/4 * \sqrt{(a + b + c) * (a - b + c) * (a + b - c) (-a + b + c)} \\
   * s &=& \frac{a + b + c}{2}
   * @f}
   */

  Real a, b, c;
  a = distance_2d(coord1, coord2);
  b = distance_2d(coord2, coord3);
  c = distance_2d(coord1, coord3);

  Real s;
  s = (a + b + c) * 0.5;

  return sqrt((s - a) * (s - b) * (s - c) / s);
}

/* -------------------------------------------------------------------------- */
inline Real Math::distance_3d(const Real * x, const Real * y) {
  return sqrt((y[0] - x[0])*(y[0] - x[0])
	      + (y[1] - x[1])*(y[1] - x[1])
	      + (y[2] - x[2])*(y[2] - x[2])
	      );
}

/* -------------------------------------------------------------------------- */
inline Real Math::tetrahedron_volume(const Real * coord) {
  Real xx[9],vol;
  const Real *x1 = coord,*x2 = coord+3,*x3 = coord+6,*x4 = coord+9;

  xx[0] = x2[0];xx[1] = x2[1];xx[2] = x2[2];
  xx[3] = x3[0];xx[4] = x3[1];xx[5] = x3[2];
  xx[6] = x4[0];xx[7] = x4[1];xx[8] = x4[2];
  vol = det3(xx);

  xx[0] = x1[0];xx[1] = x1[1];xx[2] = x1[2];
  xx[3] = x3[0];xx[4] = x3[1];xx[5] = x3[2];
  xx[6] = x4[0];xx[7] = x4[1];xx[8] = x4[2];
  vol -= det3(xx);

  xx[0] = x1[0];xx[1] = x1[1];xx[2] = x1[2];
  xx[3] = x2[0];xx[4] = x2[1];xx[5] = x2[2];
  xx[6] = x4[0];xx[7] = x4[1];xx[8] = x4[2];
  vol += det3(xx);

  xx[0] = x1[0];xx[1] = x1[1];xx[2] = x1[2];
  xx[3] = x2[0];xx[4] = x2[1];xx[5] = x2[2];
  xx[6] = x3[0];xx[7] = x3[1];xx[8] = x3[2];
  vol -= det3(xx);

  vol /= 6;

  return vol;
}

/* -------------------------------------------------------------------------- */
inline Real Math::tetrahedron_inradius(const Real * coord) {

  Real l12, l13, l14, l23, l24, l34;
  l12 = distance_3d(coord  , coord+3);
  l13 = distance_3d(coord  , coord+6);
  l14 = distance_3d(coord  , coord+9);
  l23 = distance_3d(coord+3, coord+6);
  l24 = distance_3d(coord+3, coord+9);
  l34 = distance_3d(coord+6, coord+9);

  Real s1, s2, s3, s4;

  s1 = (l12 + l23 + l13) * 0.5;
  s1 = sqrt(s1*(s1-l12)*(s1-l23)*(s1-l13));

  s2 = (l12 + l24 + l14) * 0.5;
  s2 = sqrt(s2*(s2-l12)*(s2-l24)*(s2-l14));

  s3 = (l23 + l34 + l24) * 0.5;
  s3 = sqrt(s3*(s3-l23)*(s3-l34)*(s3-l24));

  s4 = (l13 + l34 + l24) * 0.5;
  s4 = sqrt(s4*(s4-l13)*(s4-l34)*(s4-l14));

  Real volume = Math::tetrahedron_volume(coord);

  return 3*volume/(s1+s2+s3+s4);
}

/* -------------------------------------------------------------------------- */
inline void Math::barycenter(const Real * coord,
			     UInt nb_points, UInt spatial_dimension,
			     Real * barycenter) {
  memset(barycenter, 0, spatial_dimension * sizeof(Real));
  for (UInt n = 0; n < nb_points; ++n) {
    UInt offset = n * spatial_dimension;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      barycenter[i] += coord[offset + i] / (Real) nb_points;
    }
  }
}

