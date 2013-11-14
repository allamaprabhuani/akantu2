/**
 * @file   aka_types.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Thu Feb 17 17:27:24 2011
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
#include "aka_error.hh"
#include "aka_math.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */
#include <iomanip>

#ifndef __INTEL_COMPILER
#include <tr1/unordered_map>
#else
#include <map>
#endif

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_TYPES_HH__
#define __AKANTU_AKA_TYPES_HH__

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* maps                                                                       */
/* -------------------------------------------------------------------------- */
#ifndef __INTEL_COMPILER
template<class Key, class Ty>
struct unordered_map { typedef typename std::tr1::unordered_map<Key, Ty> type; };
#else
template<class Key, class Ty>
struct unordered_map { typedef typename std::map<Key, Ty> type; };
#endif

enum NormType {
  L_1 = 1,
  L_2 = 2,
  L_inf
};

template<typename T>
class Matrix;

/* ------------------------------------------------------------------------ */
/* Vector                                                                   */
/* ------------------------------------------------------------------------ */
template<typename T>
class Vector {
public:
  Vector() : n(0), values(NULL), wrapped(false) {}

  Vector(UInt n, T def = T()) : n(n), values(new T[n]), wrapped(false) {
    std::fill_n(values, n, def);
  }

  explicit Vector(T* data, UInt n) : n(n), values(data), wrapped(true) {}

  Vector(const Vector & src) {
    wrapped = src.wrapped;
    n = src.n;
    if (src.wrapped) {
      values = src.values;
    } else {
      values = new T[n];
      memcpy(this->values, src.values, n * sizeof(T));
    }
  }

  virtual ~Vector() { if(!wrapped) delete [] values; };

  /* ---------------------------------------------------------------------- */
  UInt size() const { return n; }

  T * storage() const { return values; }

  /* ---------------------------------------------------------------------- */
  void shallowCopy(const Vector & src) {
    if(!wrapped) delete [] values;
    this->n = src.n;
    this->wrapped = true;
    this->values = src.values;
  }

  /* ---------------------------------------------------------------------- */
  inline T& operator()(UInt i) { return *(values + i); };
  inline const T& operator()(UInt i) const { return *(values + i); };
  inline T& operator[](UInt idx) { return *(values + idx); };
  inline const T& operator[](UInt idx) const { return *(values + idx); };

  /* ---------------------------------------------------------------------- */
  inline Vector & operator=(const Vector & src) {
    if(this != &src) {
      if (wrapped) {
	AKANTU_DEBUG_ASSERT(n == src.n, "vectors of different size");
	memcpy(this->values, src.values, n * sizeof(T));
      } else {
	n = src.n;
	delete [] values;
	values = new T[n];
	memcpy(this->values, src.values, n * sizeof(T));
      }
    }
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator+=(const Vector & vect) {
    T * a = this->storage();
    T * b = vect.storage();
    for (UInt i = 0; i < n; ++i) *(a++) += *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator+=(const T & x) {
    T * a = this->values;
    for (UInt i = 0; i < n; ++i) *(a++) += x;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator-=(const Vector & vect) {
    T * a = this->storage();
    T * b = vect.storage();
    for (UInt i = 0; i < n; ++i) *(a++) -= *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator*=(const T & scalar) {
    T * a = this->storage();
    for (UInt i = 0; i < n; ++i) *(a++) *= scalar;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator*=(const Vector & vect) {
    T * a = this->storage();
    T * b = vect.storage();
    for (UInt i = 0; i < n; ++i) *(a++) *= *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator/=(const T & x) {
    T * a = this->values;
    for (UInt i = 0; i < n; ++i) *(a++) /= x;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Vector & operator/=(const Vector & vect) {
    T * a = this->storage();
    T * b = vect.storage();
    for (UInt i = 0; i < n; ++i) *(a++) /= *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline bool operator==(const Vector & vect) {
    T * a = this->storage();
    T * b = vect.storage();
    UInt i = 0;
    while (i < n && *(a++) == *(b++)) ++i;
    return i == n;
  }

  /* ---------------------------------------------------------------------- */
  inline Real dot(const Vector & vect) const {
    return Math::vectorDot(values, vect.storage(), n);
  }
  /* ---------------------------------------------------------------------- */
  inline Vector & crossProduct(const Vector & v1, const Vector & v2) {
    AKANTU_DEBUG_ASSERT(n == 3,
			"crossProduct is only defined in 3D (n=" << n << ")");

    AKANTU_DEBUG_ASSERT(n == v1.n && n == v2.n,
			"crossProduct is not a valid operation non matching size vectors");
    // for (UInt i = 0; i < n; ++i) {
    //	values[i] =
    //	  v1((i+1) % n) * v2((i+2) % n) -
    //	  v1((i+2) % n) * v1((i+1) % n);
    // }
    Math::vectorProduct3(v1.values, v2.values, this->values);

    return *this;
  }

  /* ------------------------------------------------------------------------ */
  inline void solve(Matrix<T> & A, const Vector<T> & b) {
    AKANTU_DEBUG_ASSERT(n == A.rows() && n = A.cols(),
			"The solution vector as a mismatch in size with the matrix");
    AKANTU_DEBUG_ASSERT(n == b.n, "The rhs vector as a mismatch in size with the matrix");
    Math::solve(n, A.storage(), values, b.storage());
  }

  /* ---------------------------------------------------------------------- */
  inline void clear() { memset(values, 0, n * sizeof(T)); };
  inline void set(T t) { std::fill_n(values, n, t); };

  template<bool tr_A>
  inline void mul(const Matrix<T> & A, const Vector & x, Real alpha = 1.0);

  /* ---------------------------------------------------------------------- */
  inline Real norm() const {
    return Math::norm(this->n, this->values);
  }

  /* ---------------------------------------------------------------------- */
  inline void normalize() {
    Real n = norm();
    operator/=(n);
  }

  /* ---------------------------------------------------------------------- */
  inline void copy(const Vector & src) {
    memcpy(values, src.storage(), n * sizeof(T));
  }

  /* ---------------------------------------------------------------------- */
  /// norm of (*this - x)
  inline Real distance(const Vector & y) const {
    Real * vx = values; Real * vy = y.storage();
    Real sum_2 = 0;
    for (UInt i = 0; i < n; ++i, ++vx, ++vy) sum_2 += (*vx - *vy)*(*vx - *vy);
    return sqrt(sum_2);
  }

  inline bool equal(Vector v, Real tolerance = Math::getTolerance()) const {
    for (UInt i = 0; i < n; ++i) {
      if(std::abs(values[i] - v(i)) > tolerance) return false;
    }
    return true;
  }


  inline short compare(Vector v, Real tolerance = Math::getTolerance()) const {
    for (UInt i(0); i < n; ++i) {
      if(std::abs(values[i] - v(i)) > tolerance)
        return values[i] - v(i) > tolerance ? 1 : -1;
    }
    return 0;
  }

  inline bool operator==(Vector v) const { return equal(v); }
  inline bool operator<(Vector v) const { return compare(v) == -1; }
  inline bool operator>(Vector v) const { return compare(v) == 1; }

  /* ---------------------------------------------------------------------- */
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << "Vector<" << debug::demangle(typeid(T).name()) << ">(" << n <<"): [";
    for (UInt i = 0; i < n; ++i) {
      if(i != 0) stream << ", ";
      stream << values[i];
    }
    stream << "]";
  }

  friend class ::akantu::Array<T>;
protected:
  UInt n;
  T * values;
  bool wrapped;
};

typedef Vector<Real> RVector;

// support operations for the creation of other vectors
template <typename T> Vector<T> operator*(T scalar, const Vector<T>& a);
template <typename T> Vector<T> operator+(const Vector<T>& a, const Vector<T>& b);
template <typename T> Vector<T> operator-(const Vector<T>& a, const Vector<T>& b);

/* -------------------------------------------------------------------------- */
template <typename T>
Vector<T> operator*(T scalar, const Vector<T>& a) {
  Vector<T> r(a.size());
  r = a;
  r *= scalar;
  return r;
}

template <typename T>
Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
  Vector<T> r(a.size());
  r = a;
  r += b;
  return r;
}

template <typename T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
  Vector<T> r(a.size());
  r = a;
  r -= b;
  return r;
}

template <typename T>
Matrix<T> operator*(T scalar, const Matrix<T>& a) {
  Matrix<T> r(a.rows(), a.cols());
  r = a;
  r *= scalar;
  return r;
}

template <typename T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {
  Matrix<T> r(a.rows(), a.cols());
  r = a;
  r += b;
  return r;
}

template <typename T>
Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b) {
  Matrix<T> r(a.rows(), a.cols());
  r = a;
  r -= b;
  return r;
}



/* ------------------------------------------------------------------------ */
/* Matrix                                                                   */
/* ------------------------------------------------------------------------ */
template<typename T>
class Matrix {
public:
  Matrix() : m(0), n(0), values(NULL), wrapped(false) {};

  Matrix(UInt m, UInt n, T def = 0) : m(m), n(n), values(new T[n*m]), wrapped(false) {
    std::fill_n(values, n*m, def);
  };

  Matrix(T * data, UInt m, UInt n) : m(m), n(n), values(data), wrapped(true) {};

  Matrix(const Matrix & src) {
    wrapped = src.wrapped;
    n = src.n;
    m = src.m;
    if (src.wrapped) {
      values = src.values;
    } else {
      values = new T[n*m];
      memcpy(this->values, src.values, n*m * sizeof(T));
    }
  }

  virtual ~Matrix() { if(!wrapped) delete [] values; };

  /* ---------------------------------------------------------------------- */
  UInt size() const { return n*m; };

  UInt rows() const { return m; };
  UInt cols() const { return n; };
  T * storage() const { return values; };

  /* ---------------------------------------------------------------------- */
  void shallowCopy(const Matrix & src) {
    if(!wrapped) delete [] values;
    this->n = src.n;
    this->m = src.m;
    this->wrapped = true;
    this->values = src.values;
  }

  /* ---------------------------------------------------------------------- */
  inline T& operator()(UInt i, UInt j) { return *(values + i + j*m); };
  inline const T& operator()(UInt i, UInt j) const { return *(values + i + j*m); };

  /// give a line vector wrapped on the column i
  inline Vector<T> operator()(UInt j) { return Vector<T>(values + j*m, m); }
  inline const Vector<T> operator()(UInt j) const { return Vector<T>(values + j*m, m); }

  inline T& operator[](UInt idx) { return *(values + idx); };
  inline const T& operator[](UInt idx) const { return *(values + idx); };
  inline Matrix & operator=(const Matrix & src) {
    if(this != &src) {
      if (wrapped) {
	AKANTU_DEBUG_ASSERT(n == src.n && m == src.m, "vectors of different size");
	memcpy(this->values, src.values, n*m * sizeof(T));
      } else {
	n = src.n;
	m = src.m;
	delete [] values;
	values = new T[n*m];
	memcpy(this->values, src.values, n*m * sizeof(T));
      }
    }
    return *this;
  };

  /* ---------------------------------------------------------------------- */
  // inline Matrix & operator=(T x) {
  //   T * a = this->values;
  //   for (UInt i = 0; i < n*m; ++i) *(a++) = x;
  //   return *this;
  // }

  /* ---------------------------------------------------------------------- */
  inline Matrix operator* (const Matrix & B) {
    Matrix C(this->m, B.n);
    C.mul<false, false>(*this, B);
    return C;
  };

  /* ----------------------------------------------------------------------- */
  inline Matrix & operator*= (const Matrix & B) {
    Matrix C(this->m, this->n);
    C.mul<false, false>(*this, B);
    this->copy(C);
    return *this;
  };

  /* ---------------------------------------------------------------------- */
  inline Matrix & operator+=(const Matrix & A) {
    T * a = this->storage();
    T * b = A.storage();
    for (UInt i = 0; i < n*m; ++i)
      *(a++) += *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix & operator-=(const Matrix & A) {
    T * a = this->storage();
    T * b = A.storage();
    for (UInt i = 0; i < n*m; ++i)
      *(a++) -= *(b++);
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix & operator+=(T x) {
    T * a = this->values;
    for (UInt i = 0; i < n*m; ++i) *(a++) += x;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix & operator*=(T x) {
    T * a = this->storage();
    for (UInt i = 0; i < n*m; ++i) *(a++) *= x;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix & operator/=(T x) {
    T * a = this->values;
    for (UInt i = 0; i < n*m; ++i) *(a++) /= x;
    return *this;
  }

  /* ---------------------------------------------------------------------- */
  template<bool tr_A, bool tr_B>
  inline void mul(const Matrix & A, const Matrix & B, T alpha = 1.0) {
    UInt k = A.n;
    if(tr_A) k = A.m;

#ifndef AKANTU_NDEBUG
    if (tr_B){
      AKANTU_DEBUG_ASSERT(k == B.n, "matrices to multiply have no fit dimensions");
      AKANTU_DEBUG_ASSERT(n == B.m, "matrices to multiply have no fit dimensions");
    }
    else {
      AKANTU_DEBUG_ASSERT(k == B.m, "matrices to multiply have no fit dimensions");
      AKANTU_DEBUG_ASSERT(n == B.n, "matrices to multiply have no fit dimensions");
    }
    if (tr_A){
      AKANTU_DEBUG_ASSERT(m == A.n, "matrices to multiply have no fit dimensions");
    }
    else{
      AKANTU_DEBUG_ASSERT(m == A.m, "matrices to multiply have no fit dimensions");
    }
#endif //AKANTU_NDEBUG

    Math::matMul<tr_A, tr_B>(m, n, k, alpha, A.storage(), B.storage(), 0., values);
  }

  /* ---------------------------------------------------------------------- */
  inline void outerProduct(const Vector<T> & A,
			   const Vector<T> & B) {
    AKANTU_DEBUG_ASSERT(A.size() == m && B.size() == n, "A and B are not compatible with the size of the matrix");
    for (UInt i = 0; i < m; ++i) {
      for (UInt j = 0; j < n; ++j) {
	values[i + j * m] += A[i] * B[j];
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  inline void eig(Vector<T> & eigenvalues, Matrix & eigenvectors) const {
    AKANTU_DEBUG_ASSERT(n == m, "eig is not a valid operation on a rectangular matrix");
    AKANTU_DEBUG_ASSERT(eigenvalues.size() == n, "eigenvalues should be of size "
                        << n << ".");
    AKANTU_DEBUG_ASSERT((eigenvectors.rows() == eigenvectors.cols()) &&
                        (eigenvectors.rows() == n),
                        "Eigenvectors needs to be a square matrix of size "
                        << n << " x " << n << ".");

    Matrix<T> temp = *this;
    Math::matrixEig(temp.n, temp.values, eigenvalues.storage(), eigenvectors.storage());
  }

  /* ---------------------------------------------------------------------- */
  inline void eig(Vector<T> & eigenvalues) const {
    AKANTU_DEBUG_ASSERT(n == m, "eig is not a valid operation on a rectangular matrix");
    AKANTU_DEBUG_ASSERT(eigenvalues.size() == n, "eigenvalues should be of size "
                        << n << ".");
    Matrix<T> temp = *this;
    Math::matrixEig(temp.n, temp.values, eigenvalues.storage());
  }

  /* ---------------------------------------------------------------------- */
  inline void clear() { std::fill_n(values, m * n, 0); };

  /* ---------------------------------------------------------------------- */
  inline void copy(const Matrix & src) {
    memcpy(values, src.storage(), m * n * sizeof(T));
  }

  /* ---------------------------------------------------------------------- */
  inline void eye(T alpha = 1.) {
    AKANTU_DEBUG_ASSERT(n == m, "eye is not a valid operation on a rectangular matrix");
    clear();
    for (UInt i = 0; i < n; ++i) {
      values[i + i * m] = alpha;
    }
  }

  /* ---------------------------------------------------------------------- */
  static inline Matrix<T> eye(UInt m, T alpha = 1.) {
    Matrix<T> tmp(m, m);
    tmp.clear();
    for (UInt i = 0; i < m; ++i) tmp(i, i) = alpha;
    return tmp;
  }

  /* ---------------------------------------------------------------------- */
  inline T trace() const {
    AKANTU_DEBUG_ASSERT(n == m, "trace is not a valid operation on a rectangular matrix");
    T trace = 0.;
    for (UInt i = 0; i < n; ++i)
      trace += values[i + i * m];
    return trace;
  }

  /* ---------------------------------------------------------------------- */
  inline Matrix transpose() const {
    Matrix tmp(n, m);
    for (UInt i = 0; i < m; ++i) {
      for (UInt j = 0; j < n; ++j) {
	tmp(j,i) = operator()(i, j);
      }
    }
    return tmp;
  }

  /* ---------------------------------------------------------------------- */
  inline void inverse(const Matrix & A) {
    AKANTU_DEBUG_ASSERT(A.n == A.m, "inv is not a valid operation on a rectangular matrix");
    AKANTU_DEBUG_ASSERT(n == A.n, "the matrix should have the same size as its inverse");

    if(n == 1) *this->values = 1./ *A.values;
    else if(n == 2) Math::inv2(A.values, this->values);
    else if(n == 3) Math::inv3(A.values, this->values);
    else Math::inv(n, A.values, this->values);
  }

  /* --------------------------------------------------------------------- */
  inline T det() {
    AKANTU_DEBUG_ASSERT(n == m, "inv is not a valid operation on a rectangular matrix");
    if(n == 1) return *values;
    else if(n == 2) return Math::det2(values);
    else if(n == 3) return Math::det3(values);
    else return Math::det(values, n);
  }

  /* --------------------------------------------------------------------- */
  inline T doubleDot(const Matrix<T> & other) {
     AKANTU_DEBUG_ASSERT(n == m, "doubleDot is not a valid operation on a rectangular matrix");
     if(n == 2) return Math::matrixDoubleDot22(values, other.values);
     else if(n == 3) return Math::matrixDoubleDot33(values, other.values);
     else AKANTU_DEBUG_ERROR("doubleDot is not defined for other spatial dimensions"
                             << " than 2 or 3.");
     return T();
  }

private:
  template<typename R, NormType norm_type>
  struct NormHelper {
    static R norm(const Matrix<R> & mat) {
      AKANTU_DEBUG_TO_IMPLEMENT();
    }
  };

  template<typename R>
  struct NormHelper<R, L_2> {
    static R norm(const Matrix<R> & mat) {
      return Math::norm(mat.size(), mat.storage());
    }
  };

  template<typename R>
  struct NormHelper<R, L_inf> {
    static R norm(const Matrix<R> & mat) {
      T _norm = 0.;
      T * it = mat.storage();
      T * end = mat.storage() + mat.size();
      for (; it < end; ++it) _norm = std::max(std::abs(*it), _norm);
      return _norm;
    }
  };

public:

  /*----------------------------------------------------------------------- */
  template<NormType norm_type>
  inline T norm() const { return NormHelper<T, norm_type>::norm(*this); }


  /* ---------------------------------------------------------------------- */
  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {
    std::string space;
    for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

    stream << space << "Matrix<" << debug::demangle(typeid(T).name())
	   << ">(" << n << "," << m <<"): " << "[";
    for (UInt i = 0; i < m; ++i) {
      if(i != 0) stream << ", ";
      stream << "[";
      for (UInt j = 0; j < n; ++j) {
	if(j != 0) stream << ", ";
	stream << operator()(i, j);
      }
      stream << "]";
    }
    stream << "]";
  };

  friend class ::akantu::Array<T>;
protected:
  UInt m;
  UInt n;
  T* values;
  bool wrapped;
};

typedef Matrix<Real> RMatrix;

/* ------------------------------------------------------------------------ */
template<typename T>
template<bool tr_A>
inline void Vector<T>::mul(const Matrix<T> & A,
			  const Vector<T> & x,
			  Real alpha) {
  UInt n = x.n;
#ifndef AKANTU_NDEBUG
  if (tr_A){
    AKANTU_DEBUG_ASSERT(n == A.rows(), "matrix and vector to multiply have no fit dimensions");
    AKANTU_DEBUG_ASSERT(this->n == A.cols(), "matrix and vector to multiply have no fit dimensions");
  } else {
    AKANTU_DEBUG_ASSERT(n == A.cols(), "matrix and vector to multiply have no fit dimensions");
    AKANTU_DEBUG_ASSERT(this->n == A.rows(), "matrix and vector to multiply have no fit dimensions");
  }
#endif
  Math::matVectMul<tr_A>(this->n, n, alpha, A.storage(), x.storage(), 0., values);
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline std::ostream & operator<<(std::ostream & stream, const Matrix<T> & _this)
{
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline std::ostream & operator<<(std::ostream & stream, const Vector<T> & _this)
{
  _this.printself(stream);
  return stream;
}

/* ------------------------------------------------------------------------ */
/* Tensor3                                                                  */
/* ------------------------------------------------------------------------ */
template<typename T>
class Tensor3 {
public:
public:
  Tensor3() : values(NULL), wrapped(true) { dim[0] = dim[1] = dim[2] = 0; };

  Tensor3(UInt m, UInt n, UInt p, T def = 0) : values(new T[n*m*p]), wrapped(false) {
    dim[0] = m; dim[1] = n; dim[2] = p;
    std::fill_n(values, n*m*p, def);
  };

  Tensor3(T * data, UInt m, UInt n, UInt p) : values(data), wrapped(true) {
    dim[0] = m; dim[1] = n; dim[2] = p;
  };

  Tensor3(const Tensor3 & src) {
    wrapped = src.wrapped;
    for (UInt i = 0; i < 3; ++i) dim[i] = src.dim[i];
    if (src.wrapped) {
      values = src.values;
    } else {
      values = new T[size()];
      memcpy(this->values, src.values, size() * sizeof(T));
    }
  }

  virtual ~Tensor3() { if(!wrapped) delete [] values; };
  /* ---------------------------------------------------------------------- */
  UInt size() const { return dim[0]*dim[1]*dim[2]; };
  UInt size(UInt i) const { return dim[i]; };
  T * storage() const { return values; };

  /* ---------------------------------------------------------------------- */
  void shallowCopy(const Tensor3 & src) {
    if(!wrapped) delete [] values;
    for (UInt i = 0; i < 3; ++i) dim[i] = src.dim[i];
    this->wrapped = true;
    this->values = src.values;
  }

  /* ---------------------------------------------------------------------- */
  inline T& operator()(UInt i, UInt j, UInt k) { return *(values + (k*dim[0] + i)*dim[1] + j); };
  inline const T& operator()(UInt i, UInt j, UInt k) const { return *(values + (k*dim[0] + i)*dim[1] + j); };
  inline Matrix<T> operator()(UInt k) { return Matrix<T>(values + k*dim[0]*dim[1], dim[0], dim[1]); }
  inline const Matrix<T> operator()(UInt k) const { return Matrix<T>(values + k*dim[0]*dim[1], dim[0], dim[1]); }
protected:
  UInt dim[3];
  T* values;
  bool wrapped;
};

__END_AKANTU__

#endif /* __AKANTU_AKA_TYPES_HH__ */
