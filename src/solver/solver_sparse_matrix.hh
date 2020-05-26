/**
 * @file   solver_sparse_matrix.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  ven mai 01 2020
 *
 * @brief A Documented file.
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
#include "aka_common.hh"
#include "dof_manager_default.hh"
#include "solver_vector.hh"
#if defined(AKANTU_PETSC)
#include "dof_manager_petsc.hh"
#include "sparse_matrix_petsc.hh"
#endif
/* -------------------------------------------------------------------------- */

#ifndef __SOLVER_SPARSE_MATRIX_H_
#define __SOLVER_SPARSE_MATRIX_H_

namespace aka {
template <class MatrixType> struct dof_manager_type {
  using type = akantu::DOFManagerDefault;
};

#if defined(AKANTU_PETSC)
template <> struct dof_manager_type<akantu::SparseMatrixPETSc> {
  using type = akantu::DOFManagerPETSc;
};
#endif

template <class MatrixType>
using dof_manager_t = typename dof_manager_type<MatrixType>::type;
} // namespace aka

namespace akantu {

/* -------------------------------------------------------------------------- */
class SolverSparseMatrix {
public:
  SolverSparseMatrix() = default;

  virtual ~SolverSparseMatrix() = default;

  virtual void applyBoundary(Real block_val = 1.) = 0;
  virtual void clearProfile() = 0;
  virtual void clear() = 0;
  virtual UInt getRelease() const = 0;
  virtual void resize() = 0;

  /// Equivalent of *gemv in blas
  virtual void matVecMul(const Array<Real> & x, Array<Real> & y,
                         Real alpha = 1., Real beta = 0.) const = 0;

  /// Equivalent of *gemv in blas
  virtual void matVecMul(const SolverVector & x, SolverVector & y,
                         Real alpha = 1., Real beta = 0.) const = 0;
};

/* -------------------------------------------------------------------------- */
template <class SparseMatrix,
          class DOFManagerType = aka::dof_manager_t<SparseMatrix>>
class SolverSparseMatrixTmpl : public SparseMatrix, public SolverSparseMatrix {
public:
  SolverSparseMatrixTmpl(DOFManagerType & dof_manager,
                         const MatrixType & matrix_type,
                         const ID & id = "solver_sparse_matrix");

  SolverSparseMatrixTmpl(SolverSparseMatrixTmpl & other,
                         const ID & id = "solver_sparse_matrix_copy");

  virtual ~SolverSparseMatrixTmpl() = default;

  void resize();

  /// modify the matrix to "remove" the blocked dof
  void applyBoundary(Real block_val = 1.) override;
  inline void clearProfile() override { SparseMatrix::clearProfile(); }
  inline void clear() override { SparseMatrix::clear(); }
  inline UInt getRelease() const override { return SparseMatrix::getRelease(); }

  void matVecMul(const Array<Real> & x, Array<Real> & y,
                 Real alpha = 1., Real beta = 0.) const override;

  /// Equivalent of *gemv in blas
  void matVecMul(const SolverVector & x, SolverVector & y, Real alpha = 1.,
                 Real beta = 0.) const override;

private:
  DOFManagerType & dof_manager;
};

} // namespace akantu
#endif // __SOLVER_SPARSE_MATRIX_H_
