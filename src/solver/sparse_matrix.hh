/**
 * @file   sparse_matrix.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  sparse matrix storage class (distributed assembled matrix)
 * This is a COO format (Coordinate List)
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SPARSE_MATRIX_HH__
#define __AKANTU_SPARSE_MATRIX_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
class TermsToAssemble;
class SolverVector;
class Communicator;
} // namespace akantu

namespace akantu {

class SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrix(const Communicator & communicator, const UInt m = 0,
               const UInt n = 0, const SizeType & size_type = SizeType::_local,
               const MatrixType & matrix_type = _unsymmetric,
               const ID & id = "sparse_matrix");

  SparseMatrix(const SparseMatrix & matrix, const ID & id = "sparse_matrix");

  virtual ~SparseMatrix();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// remove the existing profile
  virtual void clearProfile();

  /// set the matrix to 0
  virtual void clear() = 0;

  // /// add a non-zero element to the profile
  // virtual UInt add(UInt i, UInt j) = 0;

  // /// assemble a local matrix in the sparse one
  // virtual void add(UInt i, UInt j, Real value) = 0;

  /// save the profil in a file using the MatrixMarket file format
  virtual void saveProfile(const std::string & /* filename */) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// save the matrix in a file using the MatrixMarket file format
  virtual void saveMatrix(const std::string & /* filename */) const {
    AKANTU_TO_IMPLEMENT();
  };

  /// multiply the matrix by a coefficient
  virtual void mul(Real alpha) = 0;

  /// add matrices
  virtual void add(const SparseMatrix & matrix, Real alpha = 1.);

  /// copy the profile of another matrix
  virtual void copyProfile(const SparseMatrix & other) = 0;

  /// operator *=
  SparseMatrix & operator*=(Real alpha) {
    this->mul(alpha);
    return *this;
  }

protected:
  /// This is the revert of add B += \alpha * *this;
  virtual void addMeTo(SparseMatrix & B, Real alpha) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // /// return the values at potition i, j
  // virtual inline Real operator()(UInt /*i*/, UInt /*j*/) const {
  //   AKANTU_TO_IMPLEMENT();
  // }
  // /// return the values at potition i, j
  // virtual inline Real & operator()(UInt /*i*/, UInt /*j*/) {
  //   AKANTU_TO_IMPLEMENT();
  // }

  AKANTU_GET_MACRO(NbNonZero, nb_non_zero, UInt);
  UInt size() const {
    AKANTU_DEBUG_ASSERT(m == n, "The size is not defined if the matrix is not  "
                                "square, use cols and rows instead");
    return m;
  }

  UInt rows() const { return m; }

  UInt cols() const { return n; }

  AKANTU_GET_MACRO(MatrixType, matrix_type, const MatrixType &);

  virtual UInt getRelease() const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  ID id;

  /// Underlying communicators
  const Communicator & communicator;

  /// sparce matrix type
  MatrixType matrix_type;

  /// global size of the matrix
  Int m{0}, n{0};

  /// local size of the matrix
  Int m_local{0}, n_local{0};

  /// number of processors
  UInt nb_proc;

  /// number of non zero element
  UInt nb_non_zero;
};

// Array<Real> & operator*=(Array<Real> & vect, const SparseMatrix & mat);

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "sparse_matrix_inline_impl.cc"

#endif /* __AKANTU_SPARSE_MATRIX_HH__ */
