/**
 * @file   sparse_matrix_petsc.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Sat Feb 03 2018
 *
 * @brief  Implementation of PETSc matrix class
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
#include "sparse_matrix_petsc.hh"
#include "mpi_communicator_data.hh"
#include "vector_petsc.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(const Communicator & communicator, const UInt m,
                                     const UInt n, const SizeType & size_type,
                                     const MatrixType & matrix_type,
                                     const ID & id)
    : SparseMatrix(communicator, m, n, SizeType::_global, matrix_type, id) {
  AKANTU_DEBUG_IN();

  const auto & mpi_data =
      aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
  mpi_comm = mpi_data.getMPICommunicator();

  PETSc_call(MatCreate, mpi_comm, &mat);
  detail::PETScSetName(mat, id);

  m_local = m;
  n_local = n;

  switch (size_type) {
  case SizeType::_local: {
    PETSc_call(MatSetSizes, mat, m_local, n_local, PETSC_DETERMINE,
               PETSC_DETERMINE);
    break;
  }
  case SizeType::_global: {
    PETSc_call(MatSetSizes, mat, PETSC_DECIDE, PETSC_DECIDE, m, n);
    break;
  }
  }

  PETSc_call(MatSetFromOptions, mat);

  PETSc_call(MatSetUp, mat);

  PETSc_call(MatSetOption, mat, MAT_ROW_ORIENTED, PETSC_TRUE);
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  if (matrix_type == _symmetric) {
    PETSc_call(MatSetOption, mat, MAT_SYMMETRIC, PETSC_TRUE);
  }

  switch (size_type) {
  case SizeType::_local: {
    PETSc_call(MatGetSize, mat, &(this->m), &(this->n));
    break;
  }
  case SizeType::_global: {
    PETSc_call(MatGetLocalSize, mat, &m_local, &n_local);
    break;
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::SparseMatrixPETSc(const SparseMatrixPETSc & matrix,
                                     const ID & id)
    : SparseMatrix(matrix, id) {
  PETSc_call(MatDuplicate, matrix.mat, MAT_COPY_VALUES, &mat);
  detail::PETScSetName(mat, id);
}

/* -------------------------------------------------------------------------- */
SparseMatrixPETSc::~SparseMatrixPETSc() {
  AKANTU_DEBUG_IN();

  if (mat)
    PETSc_call(MatDestroy, &mat);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Method to save the nonzero pattern and the values stored at each position
 * @param filename name of the file in which the information will be stored
 */
void SparseMatrixPETSc::saveMatrix(const std::string & filename) const {
  AKANTU_DEBUG_IN();

  /// create Petsc viewer
  PetscViewer viewer;
  PETSc_call(PetscViewerASCIIOpen, mpi_comm, filename.c_str(), &viewer);
  PETSc_call(PetscViewerPushFormat, viewer, PETSC_VIEWER_ASCII_MATRIXMARKET);
  PETSc_call(MatView, mat, viewer);
  PETSc_call(PetscViewerPopFormat, viewer);
  PETSc_call(PetscViewerDestroy, &viewer);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// Equivalent of *gemv in blas
void SparseMatrixPETSc::matVecMul(const VectorPETSc & x, VectorPETSc & y,
                                  Real alpha, Real beta) const {
  // y = alpha A x + beta y
  VectorPETSc w(x, this->id + ":tmp");

  // w = A x
  if (release == 0) {
    PETSc_call(VecZeroEntries, w);
  } else {
    PETSc_call(MatMult, mat, x, w);
  }

  if (alpha != 1.) {
    // w = alpha w
    PETSc_call(VecScale, w, alpha);
  }

  // y = w + beta y
  PETSc_call(VecAYPX, y, beta, w);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::addMeToImpl(SparseMatrixPETSc & B, Real alpha) const {
  PETSc_call(MatAXPY, B.mat, alpha, mat, SAME_NONZERO_PATTERN);

  B.release++;
}

/* -------------------------------------------------------------------------- */
/**
 * Method to add another PETSc matrix to this PETSc matrix
 * @param matrix PETSc matrix to be added
 * @param alpha the factor specifying how many times the matrix should be added
 */
void SparseMatrixPETSc::addMeTo(SparseMatrix & B, Real alpha) const {
  if (aka::is_of_type<SparseMatrixPETSc>(B)) {
    auto & B_petsc = aka::as_type<SparseMatrixPETSc>(B);
    this->addMeToImpl(B_petsc, alpha);
  } else {
    AKANTU_TO_IMPLEMENT();
    //    this->addMeToTemplated<SparseMatrix>(*this, alpha);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * MatSetValues() generally caches the values. The matrix is ready to
 * use only after MatAssemblyBegin() and MatAssemblyEnd() have been
 * called. (http://www.mcs.anl.gov/petsc/)
 */
void SparseMatrixPETSc::applyModifications() {
  this->beginAssembly();
  this->endAssembly();
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::beginAssembly() {
  PETSc_call(MatAssemblyBegin, mat, MAT_FINAL_ASSEMBLY);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::endAssembly() {
  PETSc_call(MatAssemblyEnd, mat, MAT_FINAL_ASSEMBLY);
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::copyProfile(const SparseMatrix & other) {
  auto & A = aka::as_type<SparseMatrixPETSc>(other);

  MatDestroy(&mat);
  MatDuplicate(A.mat, MAT_DO_NOT_COPY_VALUES, &mat);
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::mul(Real alpha) {
  PETSc_call(MatScale, mat, alpha);
  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::clear() {
  PETSc_call(MatZeroEntries, mat);
  this->release++;
}

/* -------------------------------------------------------------------------- */
void SparseMatrixPETSc::clearProfile() {
  SparseMatrix::clearProfile();
  PETSc_call(MatResetPreallocation, mat);
  PETSc_call(MatSetOption, mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATIONS, PETSC_TRUE);
  //   PETSc_call(MatSetOption, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  clear();
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
