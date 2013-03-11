/**
 * @file   aka_blas_lapack.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
 *
 * @brief  Interface of the Fortran BLAS/LAPACK libraries
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

#ifndef __AKANTU_AKA_BLAS_LAPACK_HH__
#define __AKANTU_AKA_BLAS_LAPACK_HH__

/* -------------------------------------------------------------------------- */

#if defined(AKANTU_USE_BLAS) || defined(AKANTU_USE_LAPACK)
# include "aka_fortran_mangling.hh"
#endif //AKANTU_USE_BLAS

#ifdef AKANTU_USE_BLAS
extern "C" {
  //LEVEL 1
  double AKA_FC_GLOBAL(ddot, DDOT)(int *, double *, int *, double *, int *);
  //LEVEL 2
  int AKA_FC_GLOBAL(dgemv, DGEMV)(char *, int *, int *, double *, double *, int *,
                                  double *, int *, double *, double *, int *);
  //LEVEL 3
  int AKA_FC_GLOBAL(dgemm, DGEMM)(char *, char *, int *, int *, int *, double *,
                                  double *, int *, double *, int *, double *,
                                  double *, int *);
}

__BEGIN_AKANTU__

inline double aka_ddot(int *n, double *x, int *incx, double *y, int *incy) {
  return AKA_FC_GLOBAL(ddot, DDOT)(n, x, incx, y, incy);
}


inline int aka_dgemv(char *trans, int *m, int *n, double *
                     alpha, double *a, int *lda, double *x, int *incx,
                     double *beta, double *y, int *incy) {
  return AKA_FC_GLOBAL(dgemv, DGEMV)(trans, m, n, alpha, a, lda, x, incx,
                                     beta, y, incy);
}


inline int aka_dgemm(char *transa, char *transb,
                      int *m, int *n, int *k,
                      double *alpha, double *a, int *lda,
                      double *b, int *ldb,
                      double *beta, double *c, int *ldc) {
  return AKA_FC_GLOBAL(dgemm, DGEMM)(transa, transb, m, n, k, alpha, a, lda,
                                     b, ldb, beta, c, ldc);
}

__END_AKANTU__

#endif


#ifdef AKANTU_USE_LAPACK
extern "C" {
  // compute the eigenvalues/vectors
  void AKA_FC_GLOBAL(dgeev, DGEEV)(char* jobvl, char* jobvr, int* n, double* a,
                                   int* lda, double* wr, double* wi, double* vl, int* ldvl,
                                   double* vr, int* ldvr, double* work, int* lwork, int* info);

  // LU decomposition of a general matrix
  void AKA_FC_GLOBAL(dgetrf, DGETRF)(int* m, int *n,
                                     double* a, int* lda,
                                     int* ipiv, int* info);

  // generate inverse of a matrix given its LU decomposition
  void AKA_FC_GLOBAL(dgetri, DGETRI)(int* n, double* a, int* lda,
                                     int* ipiv, double* work, int* lwork, int* info);
}
#endif //AKANTU_USE_LAPACK

__BEGIN_AKANTU__

#ifdef AKANTU_USE_LAPACK
inline void aka_dgeev(char* jobvl, char* jobvr, int* n, double* a,
               int* lda, double* wr, double* wi, double* vl, int* ldvl,
               double* vr, int* ldvr, double* work, int* lwork, int* info) {
  AKA_FC_GLOBAL(dgeev, DGEEV)(jobvl, jobvr, n, a,
                              lda, wr, wi, vl, ldvl,
                              vr, ldvr, work, lwork, info);
}

inline void aka_dgetrf(int* m, int *n,
                double* a, int* lda,
                int* ipiv, int* info) {
  AKA_FC_GLOBAL(dgetrf, DGETRF)(m, n, a, lda, ipiv, info);
}

inline void aka_dgetri(int* n, double* a, int* lda,
                       int* ipiv, double* work, int* lwork, int* info) {
  AKA_FC_GLOBAL(dgetri, DGETRI)(n, a, lda, ipiv, work, lwork, info);
}
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused"
inline void aka_dgeev(char* jobvl, char* jobvr, int* n, double* a,
                      int* lda, double* wr, double* wi, double* vl, int* ldvl,
                      double* vr, int* ldvr, double* work, int* lwork, int* info) {
  AKANTU_DEBUG_ERROR("You have to compile with the support of LAPACK activated to use this function!");
}

inline void aka_dgetrf(int* m, int *n,
                double* a, int* lda,
                int* ipiv, int* info) {
  AKANTU_DEBUG_ERROR("You have to compile with the support of LAPACK activated to use this function!");
}

inline void aka_dgetri(int* n, double* a, int* lda,
                       int* ipiv, double* work, int* lwork, int* info) {
  AKANTU_DEBUG_ERROR("You have to compile with the support of LAPACK activated to use this function!");
}
#pragma GCC diagnostic pop
#endif

__END_AKANTU__



#endif /* __AKANTU_AKA_BLAS_LAPACK_HH__ */
