#include "dof_manager_petsc.hh"
#include "solver_vector_petsc.hh"
#include "solver_sparse_matrix.hh"
#include "sparse_matrix_petsc.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
SolverSparseMatrixTmpl<SparseMatrixPETSc, DOFManagerPETSc>::
    SolverSparseMatrixTmpl(DOFManagerPETSc & dof_manager,
                           const MatrixType & matrix_type, const ID & id)
    : SparseMatrixPETSc(dof_manager.getCommunicator(),
                        dof_manager.getPureLocalSystemSize(),
                        dof_manager.getPureLocalSystemSize(), SizeType::_local,
                        matrix_type, id),
      dof_manager(dof_manager) {
  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();
  PETSc_call(MatSetLocalToGlobalMapping, mat, is_ltog_mapping, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
template <>
SolverSparseMatrixTmpl<SparseMatrixPETSc, DOFManagerPETSc>::
    SolverSparseMatrixTmpl(SolverSparseMatrixPETSc & B, const ID & id)
    : SparseMatrixPETSc(B, id),
      dof_manager(B.dof_manager) {}

/* --------------------------------------------------------------------------*/
template <>
void SolverSparseMatrixTmpl<SparseMatrixPETSc, DOFManagerPETSc>::resize() {
  this->m_local = dof_manager.getPureLocalSystemSize();
  this->n_local = this->m_local;
  PETSc_call(MatSetSizes, mat, this->m_local, this->n_local, PETSC_DETERMINE,
             PETSC_DETERMINE);

  PETSc_call(MatGetSize, mat, &(this->m), &(this->n));

  auto & is_ltog_mapping = dof_manager.getISLocalToGlobalMapping();
  PETSc_call(MatSetLocalToGlobalMapping, mat, is_ltog_mapping, is_ltog_mapping);
}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixPETSc, DOFManagerPETSc>::applyBoundary(
    Real block_val) {
  AKANTU_DEBUG_IN();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();
  static int c = 0;

  // saveMatrix("before_blocked_" + std::to_string(c) + ".mtx");
  PETSc_call(MatZeroRowsColumnsLocal, mat, blocked_dofs.size(),
             blocked_dofs.storage(), block_val, nullptr, nullptr);

  // saveMatrix("after_blocked_" + std::to_string(c) + ".mtx");
  ++c;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixPETSc, DOFManagerPETSc>::matVecMul(
    const SolverVector & _x, SolverVector & _y, Real alpha, Real beta) const {
  auto & x = aka::as_type<SolverVectorPETSc>(_x);
  auto & y = aka::as_type<SolverVectorPETSc>(_y);
  SparseMatrixPETSc::matVecMul(x, y, alpha, beta);
}

} // namespace akantu
