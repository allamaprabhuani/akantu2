/**
 * @file   aka_math_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Jul 28 13:20:35 2010
 *
 * @brief  Implementation of the inline functions of the math toolkit
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline void Math::matrix_vector(UInt m, UInt n,
				Real * A,
				Real * x,
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
				Real * A,
				Real * B,
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
				 Real * A,
				 Real * B,
				 Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	      m, n, k,
	      1,
	      A, k,
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
				Real * A,
				Real * B,
				Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
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
				  Real * A,
				  Real * B,
				  Real * C) {
#ifdef AKANTU_USE_CBLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
	      m, n, k,
	      1,
	      A, k,
	      B, n,
	      0,
	      C, n);
#else
  memset(C, 0, m*n);
  for (UInt i = 0; i < m; ++i) {
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      UInt B_j = j * n;
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[l * m + i] * B[B_j + l];
      }
    }
  }
#endif
}
