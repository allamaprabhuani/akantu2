/**
 * @file   aka_types.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Feb 16 20:28:13 2011
 *
 * @brief  description of the "simple" types
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
#include "aka_vector.hh"

#ifdef AKANTU_USE_BLAS
# ifndef AKANTU_USE_CBLAS_MKL
#  include <cblas.h>
# else // AKANTU_USE_CBLAS_MKL
#  include <mkl_cblas.h>
# endif //AKANTU_USE_CBLAS_MKL
#endif //AKANTU_USE_BLAS

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_TYPES_HH__
#define __AKANTU_AKA_TYPES_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* matrix                                                                     */
/* -------------------------------------------------------------------------- */

template<typename T, UInt m, UInt n>
class TMatrix {
public:
  TMatrix() : values(new T[n*m]), wrapped(false) {};

  TMatrix(T def) : values(new T[n*m]), wrapped(false) {
    std::fill_n(values, n*m, def);
  };

  TMatrix(T *data) : values(data), wrapped(true) {};

  TMatrix(TMatrix & src) {
    values = src.values;
    wrapped = src.wrapped;
    src.wrapped = true;
  };

  virtual ~TMatrix() { if(!wrapped) delete [] values; };

  static UInt size() { return n*m; };


  inline T & operator()(UInt i, UInt j) { return *(values + i*m + j); };
  inline const T & operator()(UInt i, UInt j) const { return *(values + i*m + j); };
  inline T & operator[](UInt idx) { return *(values + idx); };
  inline const T & operator[](UInt idx) const { return *(values + idx); };

  inline TMatrix & operator=(const TMatrix & mat) {
    if(this != &mat) {
      if(values == NULL) { this->values = new T[n * m]; this->wrapped = false; }
      memcpy(this->values, &mat[0], size()*sizeof(T));
    }
    return *this;
  };

  inline void clear() { memset(values, 0, m * n * sizeof(T)); };

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << typeid(*this).name() << "," << n << "," << m <<"> :" << std::endl;
    for (UInt i = 0; i < m; ++i) {
      stream << space << space << "| ";
      for (UInt j = 0; j < n; ++j) {
	stream << values[i*n +j] << " ";
      }
      stream << "|" << std::endl;
    }
  };

  // T* &data() { return values; };

  friend class Vector<T>;

protected:
  T * values;
  bool wrapped;
};

/* -------------------------------------------------------------------------- */

template <UInt m, UInt n>
class RealTMatrix : public TMatrix<Real, m, n> {
public:
  RealTMatrix() : TMatrix<Real, m, n>() {};
  RealTMatrix(Real def) : TMatrix<Real, m, n>(def) {};
  RealTMatrix(Real * data) : TMatrix<Real, m, n>(data) {};
  RealTMatrix(RealTMatrix & mat) : TMatrix<Real, m, n>(mat) {};
};


template <UInt m, UInt n, UInt k>
inline RealTMatrix<m,n> operator* (const RealTMatrix<m,k> & A, const RealTMatrix<k,n> & B) {
  RealTMatrix<m,n> C;
  C.clear();
  for (UInt i = 0; i < m; ++i) {
    UInt A_i = i * k;
    UInt C_i = i * n;
    for (UInt j = 0; j < n; ++j) {
      for (UInt l = 0; l < k; ++l) {
	C[C_i + j] += A[A_i + l] * B[l * n + j];
      }
    }
  }
  return C;
}


/* -------------------------------------------------------------------------- */
/* vectors                                                                    */
/* -------------------------------------------------------------------------- */

template <typename T, UInt n>
class TVector : public TMatrix<Real, n, 1> {
public:
  TVector(Real * data) : TMatrix<Real, n, 1>(data) {};
};

template <UInt n>
class RealTVector : public TVector<Real, n> {
public:
  RealTVector(Real * data) : TVector<Real, n>(data) {};
};


/* -------------------------------------------------------------------------- */
/* common                                                                     */
/* -------------------------------------------------------------------------- */

template <typename T, UInt m, UInt n>
inline std::ostream & operator<<(std::ostream & stream, const TMatrix<T,n,m> & _this)
{
  _this.printself(stream);
  return stream;
}


/* -------------------------------------------------------------------------- */
class Matrix {
public:
  Matrix() : m(0), n(0), values(NULL), wrapped(true) {};

  Matrix(UInt m, UInt n, Real def = 0) : m(m), n(n), values(new Real[n*m]), wrapped(false) {
    std::fill_n(values, n*m, def);
  };

  Matrix(Real* data, UInt m, UInt n) : m(m), n(n), values(data), wrapped(true) {};

  Matrix(Matrix & src) {
    values = src.values;
    wrapped = src.wrapped;
    src.wrapped = true;
  };

  virtual ~Matrix() { if(!wrapped) delete [] values; };

  UInt size() { return n*m; };


  inline Real& operator()(UInt i, UInt j) { return *(values + i*m + j); };
  inline const Real& operator()(UInt i, UInt j) const { return *(values + i*m + j); };
  inline Real& operator[](UInt idx) { return *(values + idx); };
  inline const Real& operator[](UInt idx) const { return *(values + idx); };

  inline Matrix & operator=(const Matrix & mat) {
    if(this != &mat) {
      if(values != NULL) {
	memcpy(this->values, mat.values, m*n*sizeof(Real));
      } else {
	this->m = mat.m;
	this->n = mat.n;

	this->values = mat.values;
	const_cast<Matrix &>(mat).wrapped = true;
	this->wrapped = false;
      }
    }
    return *this;
  };

  inline Matrix operator* (const Matrix & B) {
    // UInt m = this->m, n = B.n, k = this->n;
    // Matrix C(m,n);
    // C.clear();

    // UInt A_i = 0, C_i = 0;
    // for (UInt i = 0; i < m; ++i, A_i += k, C_i += n)
    //   for (UInt j = 0; j < n; ++j)
    // 	for (UInt l = 0; l < k; ++l)
    // 	  C[C_i + j] += (*this)[A_i + l] * B[l * n + j];
    Matrix C(m,n);
    C.mul(*this, B);

    return C;
  };

  inline void mul(const Matrix & A, const Matrix & B) {
    UInt m = A.m, n = B.n, k = A.n;
#ifdef AKANTU_USE_BLAS
  ///  C := alpha*op(A)*op(B) + beta*C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	      m, n, k,
	      1,
	      A.values, k,
	      B.values, n,
	      0,
	      this->values, n);
#else
    this->clear();

    UInt A_i = 0, C_i = 0;
    for (UInt i = 0; i < m; ++i, A_i += k, C_i += n)
      for (UInt j = 0; j < n; ++j)
	for (UInt l = 0; l < k; ++l)
	  (*this)[C_i + j] += A[A_i + l] * B[l * n + j];
#endif
  };


  inline void clear() { memset(values, 0, m * n * sizeof(Real)); };

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << typeid(*this).name() << "," << n << "," << m <<"> :" << std::endl;
    for (UInt i = 0; i < m; ++i) {
      stream << space << space << "| ";
      for (UInt j = 0; j < n; ++j) {
	stream << values[i*n +j] << " ";
      }
      stream << "|" << std::endl;
    }
  };

  // Real* &data() { return values; };
  // bool &is_wrapped() { return wrapped; };
  friend class Vector<Real>;

protected:
  UInt m;
  UInt n;
  Real* values;
  bool wrapped;
};



inline std::ostream & operator<<(std::ostream & stream, const Matrix & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_AKA_TYPES_HH__ */
