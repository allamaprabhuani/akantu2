/**
 * @file   solver_petsc.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 13 2014
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Solver class implementation for the petsc solver
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solver_petsc.hh"
#include "mpi_communicator_data.hh"
#include "sparse_matrix_petsc.hh"
#include "vector_petsc.hh"
#include "petsc_wrapper.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SolverPETSc::SolverPETSc(const Communicator & communicator, const ID & id) {
  /// create a solver context
  const auto & mpi_data =
      aka::as_type<MPICommunicatorData>(communicator.getCommunicatorData());
  auto mpi_comm = mpi_data.getMPICommunicator();

  PETSc_call(KSPCreate, mpi_comm, &this->ksp);
  detail::PETScSetName(ksp, id);

  // If this is not called the solution vector is zeroed in the call to
  // KSPSolve().
  PETSc_call(KSPSetInitialGuessNonzero, ksp, PETSC_TRUE);
  PETSc_call(KSPSetFromOptions, ksp);
}

/* -------------------------------------------------------------------------- */
SolverPETSc::~SolverPETSc() {
  AKANTU_DEBUG_IN();

  if (ksp)
    PETSc_call(KSPDestroy, &ksp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::setOperators(SparseMatrixPETSc & A) {
  // set the matrix that defines the linear system and the matrix for
// preconditioning (here they are the same)
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
  PETSc_call(KSPSetOperators, ksp, A.getMat(), A.getMat());
#else
  PETSc_call(KSPSetOperators, ksp, A.getMat(), SAME_NONZERO_PATTERN);
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverPETSc::solve(VectorPETSc & solution, const VectorPETSc & rhs) {
  PETSc_call(KSPSolve, ksp, rhs, solution);
}

} // namespace akantu
