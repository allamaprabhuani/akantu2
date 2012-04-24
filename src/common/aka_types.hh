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

namespace types {
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

    Vector(T* data, UInt n) : n(n), values(data), wrapped(true) {}

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
	  delete []values;
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
    inline Vector & operator-=(const Vector & vect) {
      T * a = this->storage();
      T * b = vect.storage();
      for (UInt i = 0; i < n; ++i) *(a++) -= *(b++);
      return *this;
    }

    /* ---------------------------------------------------------------------- */
    inline Vector & operator*=(T & scalar) {
      T * a = this->storage();
      for (UInt i = 0; i < n; ++i) *(a++) *= scalar;
      return *this;
    }

    /* ---------------------------------------------------------------------- */
    inline Vector & operator/=(T & x) {
      T * a = this->values;
      for (UInt i = 0; i < n; ++i) *(a++) /= x;
      return *this;
    }

    /* ---------------------------------------------------------------------- */
    inline Real dot(const Vector & vect) {
      return Math::vectorDot(values, vect.storage(), n);
    }

    /* ---------------------------------------------------------------------- */
    inline Vector & crossProduct(const Vector & v1, const Vector & v2) {
      AKANTU_DEBUG_ASSERT(n == 3,
			  "crossProduct is only defined in 3D");

      AKANTU_DEBUG_ASSERT(n == v1.n && n == v2.n,
			  "crossProduct is not a valid operation non matching size vectors");
      // for (UInt i = 0; i < n; ++i) {
      // 	values[i] =
      // 	  v1((i+1) % n) * v2((i+2) % n) -
      // 	  v1((i+2) % n) * v1((i+1) % n);
      // }
      Math::vectorProduct3(v1.values, v2.values, this->values);

      return *this;
    }

    /* ---------------------------------------------------------------------- */
    inline void clear() { memset(values, 0, n * sizeof(T)); };

    template<bool tr_A>
    inline void mul(const Matrix & A, const Vector & x, Real alpha = 1.0);

    /* ---------------------------------------------------------------------- */
    inline Real norm() const {
      return Math::norm(this->n, this->values);
    }

    /* ---------------------------------------------------------------------- */
    /// norm of (*this - x)
    inline Real distance(const Vector & y) const {
      Real * vx = values; Real * vy = y.storage();
      Real sum_2 = 0;
      for (UInt i = 0; i < n; ++i, ++vx, ++vy) sum_2 += (*vx - *vy)*(*vx - *vy);
      return sqrt(sum_2);
    }

    /* ---------------------------------------------------------------------- */
    /// function to print the containt of the class
    virtual void printself(std::ostream & stream, int indent = 0) const {
      std::string space;
      for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

      stream << space << "types::Vector<" << debug::demangle(typeid(T).name()) << "> [" << n <<"] :" << std::endl;
      stream << space << AKANTU_INDENT << "| ";
      for (UInt i = 0; i < n; ++i) {
	stream << values[i] << " ";
      }
      stream << "|" << std::endl;
    }

    friend class ::akantu::Vector<T>;
  protected:
    UInt n;
    T * values;
    bool wrapped;
  };

  typedef Vector<Real> RVector;

  // support opertions for the creation of other vectors
  template <typename T> Vector<T> operator*(T scalar, const Vector<T>& a);
  template <typename T> Vector<T> operator+(const Vector<T>& a, const Vector<T>& b);
  template <typename T> Vector<T> operator-(const Vector<T>& a, const Vector<T>& b);

  /* -------------------------------------------------------------------------- */
  template <typename T>
  Vector<T> operator*(T scalar, const Vector<T>& a) {
    Vector<T> r = a;
    r *= scalar;
    return r;
  }

  template <typename T>
  Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
    Vector<T> r = a;
    r += b;
    return r;
  }

  template <typename T>
  Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
    Vector<T> r = a;
    r -= b;
    return r;
  }




  /* ------------------------------------------------------------------------ */
  /* Matrix                                                                   */
  /* ------------------------------------------------------------------------ */
  class Matrix {
  public:
    Matrix() : m(0), n(0), values(NULL), wrapped(true) {};

    Matrix(UInt m, UInt n, Real def = 0) : m(m), n(n), values(new Real[n*m]), wrapped(false) {
      std::fill_n(values, n*m, def);
    };

    Matrix(Real* data, UInt m, UInt n) : m(m), n(n), values(data), wrapped(true) {};

    Matrix(const Matrix & src) {
      m = src.m;
      n = src.n;
      values = src.values;
      wrapped = src.wrapped;
      const_cast<Matrix &>(src).wrapped = true;
    };

    virtual ~Matrix() { if(!wrapped) delete [] values; };

    /* ---------------------------------------------------------------------- */
    UInt size() const { return n*m; };

    UInt rows() const { return m; };
    UInt cols() const { return n; };
    Real * storage() const { return values; };

    /* ---------------------------------------------------------------------- */
    inline Real& operator()(UInt i, UInt j) { return *(values + i*n + j); };
    inline const Real& operator()(UInt i, UInt j) const { return *(values + i*n + j); };
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

    /* ---------------------------------------------------------------------- */
    inline Matrix operator* (const Matrix & B) {
      Matrix C(this->m, B.n);
      C.mul<true, false>(*this, B);

      return C;
    };

    /* ---------------------------------------------------------------------- */
    inline Matrix & operator*=(Real x) {
      Real * a = this->storage();
      for (UInt i = 0; i < n*m; ++i) *(a++) *= x;
      return *this;
    }

    /* ---------------------------------------------------------------------- */
    inline Matrix & operator/=(Real x) {
      Real * a = this->values;
      for (UInt i = 0; i < n*m; ++i) *(a++) /= x;
      return *this;
    }


    /* ---------------------------------------------------------------------- */
    template<bool tr_A, bool tr_B>
    inline void mul(const Matrix & A, const Matrix & B, Real alpha = 1.0) {
      UInt k = A.n;
      if(tr_A) k = A.m;

      Math::matMul<tr_A, tr_B>(m, n, k, alpha, A.storage(), B.storage(), 0., values);
    }

    /* ---------------------------------------------------------------------- */
    inline void eig(types::Vector<Real> & eigenvalues, Matrix & eigenvectors) const {
      AKANTU_DEBUG_ASSERT(n == m, "eig is not a valid operation on a rectangular matrix");
      Math::matrixEig(this->n, this->values, eigenvalues.storage(), eigenvectors.storage());
    }

    /* ---------------------------------------------------------------------- */
    inline void clear() { std::fill_n(values, m * n, 0); };

    /* ---------------------------------------------------------------------- */
    inline void copy(const Matrix & src) {
      memcpy(values, src.storage(), m * n * sizeof(Real));
    }

    /* ---------------------------------------------------------------------- */
    inline void eye(Real alpha = 1.) {
      AKANTU_DEBUG_ASSERT(n == m, "eye is not a valid operation on a rectangular matrix");
      clear();
      for (UInt i = 0; i < n; ++i) {
	values[i*n + i] = alpha;
      }
    }

    /* ---------------------------------------------------------------------- */
    inline Real trace() const {
      AKANTU_DEBUG_ASSERT(n == m, "trace is not a valid operation on a rectangular matrix");
      Real trace = 0.;
      for (UInt i = 0; i < n; ++i)
	trace += values[i*n + i];
      return trace;
    }

    /* ---------------------------------------------------------------------- */
    inline Matrix transpose() const {
      Matrix tmp(m, n);
      for (UInt i = 0; i < n; ++i) {
	for (UInt j = 0; j < m; ++j) {
	  tmp(j,i) = operator()(i,j);
	}
      }
      return tmp;
    }


    /* ---------------------------------------------------------------------- */
    /// function to print the containt of the class
    virtual void printself(std::ostream & stream, int indent = 0) const {
      std::string space;
      for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

      stream << space << "types::Matrix" << " [" << n << "," << m <<"] :" << std::endl;
      for (UInt i = 0; i < m; ++i) {
	stream << space << AKANTU_INDENT << "| ";
	for (UInt j = 0; j < n; ++j) {
	  stream << std::setw(10) << values[i*n +j] << " ";
	}
	stream << "|" << std::endl;
      }
    };

    friend class ::akantu::Vector<Real>;
  protected:
    UInt m;
    UInt n;
    Real* values;
    bool wrapped;
  };

  /* ------------------------------------------------------------------------ */
  template<typename T>
  template<bool tr_A>
  inline void Vector<T>::mul(const Matrix & A,
			     const Vector<T> & x,
			     Real alpha) {
    UInt n = x.n;
    Math::matVectMul<tr_A>(this->n, n, alpha, A.storage(), x.storage(), 0., values);
  }


  /* -------------------------------------------------------------------------- */
  inline std::ostream & operator<<(std::ostream & stream, const Matrix & _this)
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

}

__END_AKANTU__

#endif /* __AKANTU_AKA_TYPES_HH__ */
