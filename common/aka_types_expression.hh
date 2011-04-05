/**
 * @file   aka_types_expression.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sat Mar 26 22:20:09 2011
 *
 * @brief expressions  template for  types operations, inspired  from a  work of
 * Alejandro Aragón
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
#include "aka_types.hh"

#ifndef __AKANTU_AKA_TYPES_EXPRESSION_HH__
#define __AKANTU_AKA_TYPES_EXPRESSION_HH__

__BEGIN_AKANTU__

namespace types {

  /* ------------------------------------------------------------------------ */
  template <class E> struct ResultType;

  /* ------------------------------------------------------------------------ */
  template <class A>
  class Expression {
  public:
    typedef A Type;
    typedef typename A::Return   Return;
    typedef typename A::LeftVal  LeftVal;
    typedef typename A::RightVal RightVal;

    Expression(const A & a) : exp(a) { };

    Return operator() () { return exp; };
    LeftVal  left()  const { return exp.left(); };
    RightVal right() const { return exp.right(); };


  private:
    A & exp;
  };

  template <class L, class R, class Op>
  class BinaryOperation {
  public:
    typedef const L & LeftVal;
    typedef const R & RightVal;
    typedef typename ResultType<Expression<BinaryOperation<L,R,Op> > >::Type Return;
    typedef Op Operator;

    BinaryOperation(LeftVal a, RightVal b) : _left(a), _right(b) {};
    BinaryOperation(L & a, RightVal b) : _left(a), _right(b) {};

    LeftVal  left()  const { return _left; };
    RightVal right() const { return _right; };

    inline Return operator() () const {
      return Operator::apply(_left, _right);
    }

  private:
    LeftVal  _left;
    RightVal _right;
  };

  /* ------------------------------------------------------------------------ */
  template <class R, class Op>
  class UnaryOperation {
  public:
    typedef const R & RightVal;
    typedef const R & LeftVal;
    typedef typename ResultType<Expression<UnaryOperation<R,Op> > >::Type Return;
    typedef Op Operator;

    UnaryOperation(RightVal b) : _right(b) {};
    UnaryOperation(R & b) : _right(b) {};

    RightVal right() const { return _right; };

    inline Return operator() () const {
      return Operator::apply(_right);
    }

  private:
    RightVal _right;
  };

  /* ------------------------------------------------------------------------ */
  class MultiplyOp;
  class TransposeOp;

  typedef Matrix Mat;
  typedef Vector Vect;
  typedef Expression<BinaryOperation<Matrix, Matrix, MultiplyOp> > MatxMat;
  typedef Expression<BinaryOperation<Real, Matrix, MultiplyOp> > SxMat;
  typedef Expression<UnaryOperation<Matrix, TransposeOp> > TrMat;
  typedef Expression<UnaryOperation<Vector, TransposeOp> > TrVect;

  template <> struct ResultType<Expression<BinaryOperation<Mat,    Mat,   MultiplyOp> > > { typedef Matrix Type; };

  template <> struct ResultType<Expression<BinaryOperation<Mat,    TrMat, MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<TrMat,  TrMat, MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<TrMat,  Mat,   MultiplyOp> > > { typedef Matrix Type; };

  template <> struct ResultType<Expression<BinaryOperation<Mat,    SxMat, MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<SxMat,  SxMat, MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<SxMat,  Mat,   MultiplyOp> > > { typedef Matrix Type; };

  template <> struct ResultType<Expression<BinaryOperation<Vect,   TrVect,MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<TrVect, Vect,  MultiplyOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<BinaryOperation<Mat,    Vect,  MultiplyOp> > > { typedef Vector Type; };

  template <> struct ResultType<Expression<BinaryOperation<Real,   Mat,   MultiplyOp> > > { typedef Matrix Type; };

  template <> struct ResultType<Expression<UnaryOperation<Vect,  TransposeOp> > > { typedef Vector Type; };
  template <> struct ResultType<Expression<UnaryOperation<Mat,   TransposeOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<UnaryOperation<TrMat, TransposeOp> > > { typedef Matrix Type; };
  template <> struct ResultType<Expression<UnaryOperation<SxMat, TransposeOp> > > { typedef Matrix Type; };


  /* ------------------------------------------------------------------------ */
  class TransposeOp {
  public:
    TransposeOp() {};
  };

  /* ------------------------------------------------------------------------ */
  class MultiplyOp {
  public:
    MultiplyOp() {};

    template <class L, class R>
    static inline typename ResultType<Expression<BinaryOperation<L, R, MultiplyOp> > >::Type apply(const L & a, const R & b) {
      return apply(a, b);
    }


    /* ---------------------------------------------------------------------- */
    /* Blas 2                                                                 */
    /* ---------------------------------------------------------------------- */
    /// @f[ y = A * x @f]
    static inline Vector apply(const Mat & a, const Vect & x) {
      Vector y(x.size());
      Math::matVectMul<false>(a.rows(), a.cols(),
			1., a.storage(), x.storage(),
			0., y.storage());
      return y;
    }

    /// @f[ y = A^t * x @f]
    static inline Vector apply(const TrMat & at, const Vect & x) {
      const Mat & a = at.right();
      Vector y(x.size());
      Math::matVectMul<false>(a.rows(), a.cols(),
			      1., a.storage(), x.storage(),
			      0., y.storage());
      return y;
    }

    /* ---------------------------------------------------------------------- */
    /* Blas 3                                                                 */
    /* ---------------------------------------------------------------------- */
    /// @f[ C = A * x^t @f]
    static inline Mat apply(const Vect & x, const TrVect & yt) {
      const Vect & y = yt.right();
      Mat c(x.size(), y.size());
      Math::matMul<false, true>(x.size(), y.size(), 1,
				1., x.storage(), y.storage(),
				0., c.storage());
      return c;
    }

    /// @f[ y = A^t * x^t @f]
    static inline Mat apply(const TrVect & xt, const Vect & y) {
      const Vect & x = xt.right();
      Mat c(y.size(), x.size());
      Math::matMul<true, false>(y.size(), x.size(), 1,
				1., x.storage(), y.storage(),
				0., c.storage());
      return c;
    }


    /// @f[ C = A * B @f]
    static inline Mat apply(const Mat & a, const Mat & b) {
      Mat c(a.rows(), b.cols());
      Math::matMul<false, false>(a.rows(), b.cols(), a.cols(),
				 1., a.storage(), b.storage(),
				 0., c.storage());
      return c;
    }

    /// @f[ C = A^t * B @f]
    static inline Mat apply(const TrMat & at, const Mat & b) {
      const Mat & a = at.right();
      Mat c(a.rows(), b.cols());
      Math::matMul<true, false>(a.rows(), b.cols(), a.cols(),
				1., a.storage(), b.storage(),
				0., c.storage());
      return c;
    }

    /// @f[ C = A * B^t @f]
    static inline Mat apply(const Mat & a, const TrMat & bt) {
      const Mat & b = bt.right();
      Mat c(a.rows(), b.cols());
      Math::matMul<false, true>(a.rows(), b.cols(), a.cols(),
				1., a.storage(), b.storage(),
				0., c.storage());
      return c;
    }

    /// @f[ C = A^t * B^t @f]
    static inline Mat apply(const TrMat & at, const TrMat & bt) {
      const Mat & a = at.right();
      const Mat & b = bt.right();
      Mat c(a.rows(), b.cols());
      Math::matMul<true, true>(a.rows(), b.cols(), a.cols(),
			       1., a.storage(), b.storage(),
			       0., c.storage());
      return c;
    }
  };

  template<class A, class B>
  inline typename ResultType<Expression<BinaryOperation<A, B, MultiplyOp> > >::Type
  operator* (const A & a, const A & b) {
    return Expression<BinaryOperation<A,B,MultiplyOp> >(BinaryOperation<A,B,MultiplyOp>(a,b));
  }

} // namespace types

__END_AKANTU__

#endif /* __AKANTU_AKA_TYPES_EXPRESSION_HH__ */
