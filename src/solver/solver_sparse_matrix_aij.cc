/* -------------------------------------------------------------------------- */
#include "solver_sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include "dof_synchronizer.hh"
#include "solver_vector_default.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::
    SolverSparseMatrixTmpl(DOFManagerDefault & dof_manager,
                           const MatrixType & matrix_type, const ID & id)
    : SparseMatrixAIJ(dof_manager.getCommunicator(),
                      dof_manager.getSystemSize(),
                      dof_manager.getSystemSize(), SizeType::_global,
                      matrix_type, id),
      dof_manager(dof_manager) {}

/* -------------------------------------------------------------------------- */
template <>
SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::
    SolverSparseMatrixTmpl(SolverSparseMatrixAIJ & B, const ID & id)
    : SparseMatrixAIJ(B, id), dof_manager(B.dof_manager) {}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::resize() {
  this->m = this->n = this->dof_manager.getSystemSize();
  this->m_local = this->n_local = this->dof_manager.getPureLocalSystemSize();
}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::applyBoundary(
    Real block_val) {
  AKANTU_DEBUG_IN();

  const auto & blocked_dofs = this->dof_manager.getGlobalBlockedDOFs();
  auto begin = blocked_dofs.begin();
  auto end = blocked_dofs.end();

  auto is_blocked = [&](auto && i) -> bool {
    auto il = this->dof_manager.globalToLocalEquationNumber(i);
    return std::binary_search(begin, end, il);
  };

  for (auto && ij_a : zip(this->irn, this->jcn, a)) {
    UInt ni = std::get<0>(ij_a) - 1;
    UInt nj = std::get<1>(ij_a) - 1;

    if (is_blocked(ni) or is_blocked(nj)) {

      std::get<2>(ij_a) =
          std::get<0>(ij_a) != std::get<1>(ij_a)
              ? 0.
              : this->dof_manager.isLocalOrMasterDOF(
                    this->dof_manager.globalToLocalEquationNumber(ni))
                    ? block_val
                    : 0.;
    }
  }

  this->value_release++;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::matVecMul(
    const Array<Real> & x, Array<Real> & y, Real alpha, Real beta) const {
  auto && id = [&](auto && i) {
    return this->dof_manager.globalToLocalEquationNumber(i);
  };

  SparseMatrixAIJ::matVecMulLocal(x, y, id, id, alpha, beta);

  if (this->dof_manager.hasSynchronizer())
    this->dof_manager.getSynchronizer().reduceSynchronizeArray<AddOperation>(y);
}

/* -------------------------------------------------------------------------- */
template <>
void SolverSparseMatrixTmpl<SparseMatrixAIJ, DOFManagerDefault>::matVecMul(
    const SolverVector & _x, SolverVector & _y, Real alpha, Real beta) const {
  AKANTU_DEBUG_IN();

  auto && x = aka::as_type<SolverVectorDefault>(_x);
  auto && y = aka::as_type<SolverVectorDefault>(_y);

  auto && id = [&](auto && i) {
    return this->dof_manager.globalToLocalEquationNumber(i);
  };

  SparseMatrixAIJ::matVecMulLocal(x, y, id, id, alpha, beta);

  if (this->dof_manager.hasSynchronizer())
    this->dof_manager.getSynchronizer().reduceSynchronizeArray<AddOperation>(y);
}

} // namespace akantu
