/**
 * @file   sparse_matrix_petsc.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 06 2018
 *
 * @brief  Interface for PETSc matrices
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

#ifndef __AKANTU_PETSC_MATRIX_HH__
#define __AKANTU_PETSC_MATRIX_HH__

/* -------------------------------------------------------------------------- */
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <petscmat.h>
/* -------------------------------------------------------------------------- */

namespace akantu {
class VectorPETSc;
}

namespace akantu {

class SparseMatrixPETSc : public SparseMatrix {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SparseMatrixPETSc(const Communicator & communicator, const UInt m = 0, const UInt n = 0,
                    const SizeType & size_type = SizeType::_local,
                    const MatrixType & matrix_type = _unsymmetric,
                    const ID & id = "sparse_matrix_petsc");

  SparseMatrixPETSc(const SparseMatrixPETSc & matrix,
                    const ID & id = "sparse_matrix_petsc");

  virtual ~SparseMatrixPETSc();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set the matrix to 0
  void clear();
  void clearProfile() override;

  /// add a non-zero element to the profile
  inline UInt add(UInt i, UInt j);

  /// assemble a local matrix in the sparse one
  inline void add(UInt i, UInt j, Real value);

  inline void addLocal(UInt i, UInt j);
  inline void addLocal(UInt i, UInt j, Real val);

  template<class Rows, class Cols, class Values>
  void addLocal(const Rows & rows, const Cols & cols,
                const Values & vals);

  /// add a block of values
  template<class Rows, class Cols, class Values>
  void addValues(const Rows & is, const Cols & js,
                 const Values & values, MatrixType values_type);

  /// save the profil in a file using the MatrixMarket file format
  // void saveProfile(__attribute__((unused)) const std::string &) const
  // override {
  //   AKANTU_DEBUG_TO_IMPLEMENT();
  // }

  /// save the matrix in a file using the MatrixMarket file format
  void saveMatrix(const std::string & filename) const;

  /// multiply the matrix by a coefficient
  void mul(Real alpha);

  /// Equivalent of *gemv in blas
  void matVecMul(const VectorPETSc & x, VectorPETSc & y, Real alpha = 1.,
                 Real beta = 0.) const;

  void matVecMul(const Vec & x, Vec & y, Real alpha = 1.,
                 Real beta = 0.) const;

  /// copy the profile of a matrix
  void copyProfile(const SparseMatrix & other);

  void applyModifications();

  // void resize();

protected:
  /// This is the revert of add B += \alpha * *this;
  void addMeTo(SparseMatrix & B, Real alpha) const override;

  /// This is the specific implementation
  void addMeToImpl(SparseMatrixPETSc & B, Real alpha) const;

  void beginAssembly();
  void endAssembly();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the values at potition i, j
  virtual inline Real operator()(__attribute__((unused)) UInt i,
                                 __attribute__((unused)) UInt j) const {
    AKANTU_TO_IMPLEMENT();
  }
  /// return the values at potition i, j
  virtual inline Real & operator()(__attribute__((unused)) UInt i,
                                   __attribute__((unused)) UInt j) {
    AKANTU_TO_IMPLEMENT();
  }

  UInt getRelease() const { return release; };

  operator Mat &() { return mat; }
  operator const Mat &() const { return mat; }
  AKANTU_GET_MACRO(Mat, mat, const Mat &);
  AKANTU_GET_MACRO_NOT_CONST(Mat, mat, Mat &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// store the PETSc matrix
  Mat mat;

  /// matrix release
  UInt release{0};

  /// communicator used
  MPI_Comm mpi_comm;
};

} // namespace akantu

#include "sparse_matrix_petsc_tmpl.hh"

#endif /* __AKANTU_PETSC_MATRIX_HH__ */
