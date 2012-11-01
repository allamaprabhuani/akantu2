/**
 * @file   aka_math_tmpl.hh
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Mathilde Radiguet <mathilde.radiguet@epfl.ch>
 * @author Alejandro Marcos Aragon <alejandro.aragon@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Anciaux Guillaume <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @date   Wed Jul 28 13:20:35 2010
 *
 * @brief  Implementation of the inline functions of the math toolkit
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


__END_AKANTU__

#include <cmath>
#include <cstring>

#ifdef AKANTU_USE_BLAS
# ifndef AKANTU_USE_BLAS_MKL
#  include <cblas.h>
# else // AKANTU_USE_BLAS_MKL
#  include <mkl_cblas.h>
# endif //AKANTU_USE_BLAS_MKL
#endif //AKANTU_USE_BLAS

#ifdef AKANTU_USE_LAPACK
extern "C" {
  void dgeev_(char* jobvl, char* jobvr, int* n, double* a,
	      int* lda, double* wr, double* wi, double* vl, int* ldvl,
	      double* vr, int* ldvr, double* work, int* lwork, int* info);

  // LU decomposition of a general matrix
  void dgetrf_(int* m, int *n,
	       double* a, int* lda,
	       int* ipiv, int* info);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* n, double* a, int* lda,
	       int* ipiv, double* work, int* lwork, int* info);
}
#endif

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
inline void Math::matrix_vector(UInt m, UInt n,
				const Real * A,
				const Real * x,
				Real * y, Real alpha) {
#ifdef AKANTU_USE_BLAS
  /// y = alpha*op(A)*x + beta*y
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      m, n, alpha, A, n, x, 1, 0, y, 1);
#else
  memset(y, 0, m*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    UInt A_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      y[i] += A[A_i + j] * x[j];
    }
    y[i] *= alpha;
  }
#endif
}


/* -------------------------------------------------------------------------- */
inline void Math::matrixt_vector(UInt m, UInt n,
				 const Real * A,
				 const Real * x,
				 Real * y, Real alpha) {
#ifdef AKANTU_USE_BLAS
  /// y = alpha*op(A)*x + beta*y
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      m, n, alpha, A, m, x, 1, 0, y, 1);
#else
  memset(y, 0, m*sizeof(Real));
  for (UInt i = 0; i < m; ++i) {
    for (UInt j = 0; j < n; ++j) {
      y[i] += A[i + j * m] * x[j];
    }
    y[i] *= alpha;
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrix(UInt m, UInt n, UInt k,
				const Real * A,
				const Real * B,
				Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      m, n, k,
	      alpha,
	      A, k,
	      B, n,
	      0,
	      C, n);
#else
  memset(C, 0, m*n*sizeof(Real));
  for (UInt j = 0; j < n; ++j) {
    for (UInt i = 0; i < m; ++i) {
      UInt A_i = i * k;
      UInt C_i = i * n;
      for (UInt l = 0; l < k; ++l) {
	UInt B_l = l * n;
	C[C_i + j] += A[A_i + l] * B[B_l + j];
      }
      C[C_i + j] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrix(UInt m, UInt n, UInt k,
				 const Real * A,
				 const Real * B,
				 Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	      m, n, k,
	      alpha,
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
      C[C_i + j] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix_matrixt(UInt m, UInt n, UInt k,
				const Real * A,
				const Real * B,
				Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
	      m, n, k,
	      alpha,
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
      C[C_i + j] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline void Math::matrixt_matrixt(UInt m, UInt n, UInt k,
				  const Real * A,
				  const Real * B,
				  Real * C, Real alpha) {
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	      m, n, k,
	      alpha,
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
      C[C_i + j] *= alpha;
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
inline Real Math::vectorDot(const Real * v1, const Real * v2, UInt n) {
#ifdef AKANTU_USE_BLAS
  ///  d := v1 . v2
  Real d = cblas_ddot(n, v1, 1, v2, 1);
#else
  Real d = 0;
  for (UInt i = 0; i < n; ++i) {
    d += v1[i] * v2[i];
  }
#endif
  return d;
}


/* -------------------------------------------------------------------------- */
template <bool tr_A, bool tr_B>
inline void Math::matMul(UInt m, UInt n, UInt k,
			 Real alpha, const Real * A, const Real * B,
			 __attribute__ ((unused)) Real beta,  Real * C) {
  if(tr_A) {
    if(tr_B) matrixt_matrixt(m, n, k, A, B, C, alpha);
    else matrixt_matrix(m, n, k, A, B, C, alpha);
  } else {
    if(tr_B) matrix_matrixt(m, n, k, A, B, C, alpha);
    else matrix_matrix(m, n, k, A, B, C, alpha);
  }
}

/* -------------------------------------------------------------------------- */
template <bool tr_A>
inline void Math::matVectMul(UInt m, UInt n,
			     Real alpha, const Real * A, const Real * x,
			     __attribute__ ((unused)) Real beta, Real * y) {
  if(tr_A) {
    matrixt_vector(m, n, A, x, y, alpha);
  } else {
    matrix_vector(m, n, A, x, y, alpha);
  }
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_LAPACK
inline void Math::matrixEig(UInt n, Real * A, Real * d, Real * V) {

  // Matrix  A is  row major,  so the  lapack function  in fortran  will process
  // A^t. Asking for the left eigenvectors of A^t will give the transposed right
  // eigenvectors of A so in the C++ code the right eigenvectors.
  char jobvl;
  if(V != NULL)
    jobvl = 'V'; // compute left  eigenvectors
  else
    jobvl = 'N'; // compute left  eigenvectors

  char jobvr('N'); // compute right eigenvectors

  double * di = new double[n]; // imaginary part of the eigenvalues

  int info;
  int N = n;

  double wkopt;
  int lwork = -1;
  // query and allocate the optimal workspace
  dgeev_(&jobvl, &jobvr, &N, A, &N, d, di, V, &N, NULL, &N, &wkopt, &lwork, &info);

  lwork = int(wkopt);
  double * work = new double[lwork];
  // solve the eigenproblem
  dgeev_(&jobvl, &jobvr, &N, A, &N, d, di, V, &N, NULL, &N, work, &lwork, &info);

  AKANTU_DEBUG_ASSERT(info == 0, "Problem computing eigenvalues/vectors. DGEEV exited with the value " << info);


  delete [] work;
  delete [] di; // I hope for you that there was no complex eigenvalues !!!
}
#else
inline void Math::matrixEig(__attribute__((unused)) UInt n,
			    __attribute__((unused)) Real * A,
			    __attribute__((unused)) Real * d,
			    __attribute__((unused)) Real * V) {
  AKANTU_DEBUG_ERROR("You have to compile with the support of LAPACK activated to use this function!");
}
#endif

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_LAPACK
inline void Math::inv(UInt n, const Real * A, Real * invA) {
  int N = n;
  int info;
  int * ipiv = new int[N+1];
  int lwork = N*N;
  double * work = new double[lwork];

  std::copy(A, A + n*n, invA);

  dgetrf_(&N, &N, invA, &N, ipiv, &info);
  if(info > 0) {
    AKANTU_DEBUG_ERROR("Singular matrix - cannot factorize it (info: "
		       << info <<" )");
  }

  dgetri_(&N, invA, &N, ipiv, work, &lwork, &info);
  if(info != 0) {
    AKANTU_DEBUG_ERROR("Cannot invert the matrix (info: "<< info <<" )");
  }

  delete [] ipiv;
  delete [] work;
}
#else
inline void Math::inv(__attribute__((unused)) UInt n,
		      __attribute__((unused)) const Real * A,
		      __attribute__((unused)) Real * Ainv) {
  AKANTU_DEBUG_ERROR("You have to compile with the support of LAPACK activated to use this function!");
}
#endif

/* -------------------------------------------------------------------------- */
inline void Math::matrix22_eigenvalues(Real * A, Real *Adiag) {
  ///d = determinant of Matrix A
  Real d = det2(A);
  ///b = trace of Matrix A
  Real b = A[0]+A[3];

  Real c = sqrt(b*b - 4 *d);
  Adiag[0]= .5*(b + c);
  Adiag[1]= .5*(b - c);
}

/* -------------------------------------------------------------------------- */
inline void Math::matrix33_eigenvalues(Real * A, Real *Adiag) {
  /// a L^3 + b L^2 + c L + d = 0
  matrixEig(3, A, Adiag);
//  Real a = -1 ;
//  ///b = trace of Matrix A
//  Real b = A[0]+A[4]+A[8];
//  /// c = 0.5*(trace(M^2)-trace(M)^2)
//  Real c =  A[1]*A[3] + A[2]*A[6] + A[5]*A[7] - A[0]*A[4] -
//    A[0]*A[8] - A[4]*A[8];
//  ///d = determinant of Matrix A
//  Real d = det3(A);
//
//  /// Define x, y, z
//  Real x = c/a - b*b/(3.*a*a);
//  Real y = 2.*b*b*b/(27.*a*a*a) - b*c/(3.*a*a) + d/a;
//  Real z = y*y/4. + x*x*x/27.;
//  /// Define I, j, k, m, n, p (so equations are not so cluttered)
//  Real i = sqrt(y*y/4. - z);
//  Real j = pow(i,1./3.);
//  Real k = 0;
//  if (std::abs(i) > 1e-12)
//    k = acos(-(y/(2.*i)));
//
//  Real m = cos(k/3);
//  Real n = sqrt(3.)*sin(k/3);
//  Real p = -b/(3.*a);
//
//  Adiag[0] = 2*j*m + p;
//  Adiag[1] = -j *(m + n) + p;
//  Adiag[2] = -j * (m - n) + p;
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Math::eigenvalues(Real * A, Real * d) {
  if(dim == 1) { d[0] = A[0]; }
  else if(dim == 2) { matrix22_eigenvalues(A, d); }
  // else if(dim == 3) { matrix33_eigenvalues(A, d); }
  else matrixEig(dim, A, d);
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
inline Real Math::norm(UInt n, const Real * vec) {
  Real norm = 0.;
  for (UInt i = 0; i < n; ++i) {
    norm += vec[i]*vec[i];
  }
  return sqrt(norm);
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
inline Real Math::vectorDot2(const Real * v1, const Real * v2) {
  return (v1[0]*v2[0] + v1[1]*v2[1]);
}

/* -------------------------------------------------------------------------- */
inline Real Math::vectorDot3(const Real * v1, const Real * v2) {
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
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
inline Real Math::tetrahedron_volume(const Real * coord1,
				     const Real * coord2,
				     const Real * coord3,
				     const Real * coord4) {
  Real xx[9], vol;

  xx[0] = coord2[0]; xx[1] = coord2[1]; xx[2] = coord2[2];
  xx[3] = coord3[0]; xx[4] = coord3[1]; xx[5] = coord3[2];
  xx[6] = coord4[0]; xx[7] = coord4[1]; xx[8] = coord4[2];
  vol = det3(xx);

  xx[0] = coord1[0]; xx[1] = coord1[1]; xx[2] = coord1[2];
  xx[3] = coord3[0]; xx[4] = coord3[1]; xx[5] = coord3[2];
  xx[6] = coord4[0]; xx[7] = coord4[1]; xx[8] = coord4[2];
  vol -= det3(xx);

  xx[0] = coord1[0]; xx[1] = coord1[1]; xx[2] = coord1[2];
  xx[3] = coord2[0]; xx[4] = coord2[1]; xx[5] = coord2[2];
  xx[6] = coord4[0]; xx[7] = coord4[1]; xx[8] = coord4[2];
  vol += det3(xx);

  xx[0] = coord1[0]; xx[1] = coord1[1]; xx[2] = coord1[2];
  xx[3] = coord2[0]; xx[4] = coord2[1]; xx[5] = coord2[2];
  xx[6] = coord3[0]; xx[7] = coord3[1]; xx[8] = coord3[2];
  vol -= det3(xx);

  vol /= 6;

  return vol;
}

/* -------------------------------------------------------------------------- */
inline Real Math::tetrahedron_inradius(const Real * coord1,
				       const Real * coord2,
				       const Real * coord3,
				       const Real * coord4) {

  Real l12, l13, l14, l23, l24, l34;
  l12 = distance_3d(coord1, coord2);
  l13 = distance_3d(coord1, coord3);
  l14 = distance_3d(coord1, coord4);
  l23 = distance_3d(coord2, coord3);
  l24 = distance_3d(coord2, coord4);
  l34 = distance_3d(coord3, coord4);

  Real s1, s2, s3, s4;

  s1 = (l12 + l23 + l13) * 0.5;
  s1 = sqrt(s1*(s1-l12)*(s1-l23)*(s1-l13));

  s2 = (l12 + l24 + l14) * 0.5;
  s2 = sqrt(s2*(s2-l12)*(s2-l24)*(s2-l14));

  s3 = (l23 + l34 + l24) * 0.5;
  s3 = sqrt(s3*(s3-l23)*(s3-l34)*(s3-l24));

  s4 = (l13 + l34 + l14) * 0.5;
  s4 = sqrt(s4*(s4-l13)*(s4-l34)*(s4-l14));

  Real volume = Math::tetrahedron_volume(coord1,coord2,coord3,coord4);

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

/* -------------------------------------------------------------------------- */
inline void Math::vector_2d(const Real * x, const Real * y, Real * res) {
  res[0] = y[0]-x[0];
  res[1] = y[1]-x[1];
}

/* -------------------------------------------------------------------------- */
inline void Math::vector_3d(const Real * x, const Real * y, Real * res) {
  res[0] = y[0]-x[0];
  res[1] = y[1]-x[1];
  res[2] = y[2]-x[2];
}

/* -------------------------------------------------------------------------- */
inline bool Math::are_float_equal(const Real x, const Real y){
  return (std::abs( x - y) < tolerance);
}
/* -------------------------------------------------------------------------- */
inline bool Math::isnan(Real x) {
#if defined(__INTEL_COMPILER)
#pragma warning ( push )
#pragma warning ( disable : 1572 )
#endif //defined(__INTEL_COMPILER)

  // x = x return false means x = quiet_NaN
  return !(x == x); 

#if defined(__INTEL_COMPILER)
#pragma warning ( pop )
#endif //defined(__INTEL_COMPILER)
}

/* -------------------------------------------------------------------------- */
inline bool Math::are_vector_equal(UInt n, Real * x, Real * y){
  bool test = true;
  for (UInt i = 0; i < n; ++i) {
    test &= are_float_equal(x[i],y[i]);
  }

  return test;
}

/* -------------------------------------------------------------------------- */
inline bool Math::intersects(Real x_min, Real x_max, Real y_min, Real y_max) {
  return ! ((x_max < y_min) || (x_min > y_max));
}

/* -------------------------------------------------------------------------- */
inline bool Math::is_in_range(Real a, Real x_min, Real x_max) {
  return ((a >= x_min) && (a <= x_max));
}
