/* -------------------------------------------------------------------------- */
#include "sparse_matrix_petsc.hh"
#include "petsc_wrapper.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt SparseMatrixPETSc::add(UInt i, UInt j) {
  PETSc_call(MatSetValue, mat, i, j, 0, ADD_VALUES);
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixPETSc::add(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValue, mat, i, j, val, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixPETSc::addLocal(UInt i, UInt j) {
  PETSc_call(MatSetValueLocal, mat, i, j, 0, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
inline void SparseMatrixPETSc::addLocal(UInt i, UInt j, Real val) {
  PETSc_call(MatSetValueLocal, mat, i, j, val, ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
template <class Rows, class Cols, class Values>
void SparseMatrixPETSc::addLocal(const Rows & rows, const Cols & cols,
                                 const Values & vals) {
  PETSc_call(MatSetValuesLocal, mat, rows.size(), rows.storage(), cols.size(),
             cols.storage(), vals.storage(), ADD_VALUES);
}

/* -------------------------------------------------------------------------- */
template <class Rows, class Cols, class Values>
void SparseMatrixPETSc::addValues(const Rows & rows, const Cols & cols,
                                  const Values & vals, MatrixType type) {
  if (type == _unsymmetric and matrix_type == _symmetric) {
    PETSc_call(MatSetOption, mat, MAT_SYMMETRIC, PETSC_FALSE);
    PETSc_call(MatSetOption, mat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
  }

  PETSc_call(MatSetValues, mat, rows.size(), rows.storage(), cols.size(),
             cols.storage(), vals.storage(), ADD_VALUES);
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
