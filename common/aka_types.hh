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
#include "aka_math.hh"
#include "aka_vector.hh"

/* -------------------------------------------------------------------------- */
#ifndef __INTEL_COMPILER
#include <tr1/unordered_map>
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
      C.mul<false, false>(*this, B);

      return C;
    };

    template<bool tr_A, bool tr_B>
    inline void mul(const Matrix & A, const Matrix & B, Real alpha = 1.0) {
      UInt k = A.n;
      if(tr_A) k = A.m;

      Math::matMul<tr_A, tr_B>(m, n, k, alpha, A.storage(), B.storage(), 0., values);
    }

    inline void clear() { memset(values, 0, m * n * sizeof(Real)); };

    /// function to print the containt of the class
    virtual void printself(std::ostream & stream, int indent = 0) const {
      std::string space;
      for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

      stream << space << debug::demangle(typeid(*this).name()) << "<" << n << "," << m <<"> :" << std::endl;
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

  class Tensor2 : public Matrix {
  };

  /* ------------------------------------------------------------------------ */
  /* Vector                                                                   */
  /* ------------------------------------------------------------------------ */
  template<typename T>
  class Vector {
  public:
    Vector() : n(0), values(NULL), wrapped(true) {};

    Vector(UInt n, Real def = 0) : n(n), values(new T[n]), wrapped(false) {
      std::fill_n(values, n, def);
    };

    Vector(T* data, UInt n) : n(n), values(data), wrapped(true) {};

    Vector(const Vector & src) {
      n = src.n;
      values = src.values;
      wrapped = src.wrapped;
      const_cast<Vector &>(src).wrapped = true;
    };

    virtual ~Vector() { if(!wrapped) delete [] values; };

    /* ---------------------------------------------------------------------- */
    UInt size() const { return n; };

    T * storage() const { return values; };

    /* ---------------------------------------------------------------------- */
    inline T& operator()(UInt i) { return *(values + i); };
    inline const T& operator()(UInt i) const { return *(values + i); };
    inline T& operator[](UInt idx) { return *(values + idx); };
    inline const T& operator[](UInt idx) const { return *(values + idx); };

    /* ---------------------------------------------------------------------- */
    inline Vector & operator=(const Vector & vect) {
      if(this != &vect) {
	if(values != NULL) {
	  memcpy(this->values, vect.values, n * sizeof(T));
	} else {
	  this->n = vect.n;

	  this->values = vect.values;
	  const_cast<Vector &>(vect).wrapped = true;
	  this->wrapped = false;
	}
      }
      return *this;
    };

    inline void clear() { memset(values, 0, n * sizeof(T)); };

    template<bool tr_A>
    inline void mul(const Matrix & A, const Vector & x, Real alpha = 1.0) {
      UInt n = x.n;
      Math::matVectMul<tr_A>(this->n, n, alpha, A.storage(), x.storage(), 0., values);
    }

    /* ---------------------------------------------------------------------- */
    /// function to print the containt of the class
    virtual void printself(std::ostream & stream, int indent = 0) const {
      std::string space;
      for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

      stream << space << debug::demangle(typeid(*this).name()) << "<" << n <<"> :" << std::endl;
      stream << space << AKANTU_INDENT << "| ";
      for (UInt i = 0; i < n; ++i) {
	stream << values[i] << " ";
      }
      stream << "|" << std::endl;
    };

    friend class ::akantu::Vector<T>;
  protected:
    UInt n;
    T * values;
    bool wrapped;
  };

  typedef Vector<Real> RVector;


  /* -------------------------------------------------------------------------- */
  inline std::ostream & operator<<(std::ostream & stream, const Matrix & _this)
  {
    _this.printself(stream);
    return stream;
  }

  /* -------------------------------------------------------------------------- */
  inline std::ostream & operator<<(std::ostream & stream, const RVector & _this)
  {
    _this.printself(stream);
    return stream;
  }


  /* -------------------------------------------------------------------------- */
  /* templated matrix                                                           */
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


    inline T & operator()(UInt i, UInt j) { return *(values + i*n + j); };
    inline const T & operator()(UInt i, UInt j) const { return *(values + i*n + j); };
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

      stream << space << debug::demangle(typeid(*this).name())<< ":" << std::endl;
      for (UInt i = 0; i < m; ++i) {
	stream << space << AKANTU_INDENT << "| ";
	for (UInt j = 0; j < n; ++j) {
	  stream << std::setw(10) << values[i*n +j] << " ";
	}
	stream << "|" << std::endl;
      }
    };

    // T* &data() { return values; };
    friend class ::akantu::Vector<T>;

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
  /* templated vectors                                                          */
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
}

__END_AKANTU__

#endif /* __AKANTU_AKA_TYPES_HH__ */
