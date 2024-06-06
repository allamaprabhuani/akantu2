/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "sparse_solver_petsc.hh"
#include "dof_manager_petsc.hh"
#include "mpi_communicator_data.hh"
#include "solver_vector_petsc.hh"
#include "sparse_matrix_petsc.hh"
/* -------------------------------------------------------------------------- */
#include <petscksp.h>
// #include <petscsys.h>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
SparseSolverPETSc::SparseSolverPETSc(DOFManager & dof_manager,
                                     const ID & matrix_id, const ID & id)
    : SparseSolver(dof_manager, matrix_id, id),
      matrix(getDOFManager().getMatrix(matrix_id)) {
  auto && mpi_comm = getDOFManager().getMPIComm();

  this->registerParam("petsc_options", petsc_options, _pat_parsable,
                      "PETSc options");

  /// create a solver context
  PETSc_call(KSPCreate, mpi_comm, &this->ksp);
}
/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::initialize() {}
/* -------------------------------------------------------------------------- */
SparseSolverPETSc::~SparseSolverPETSc() {
  if (ksp != nullptr) {
    PETSc_call(KSPDestroy, &ksp);
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::setOperators() {
  // set the matrix that defines the linear system and the matrix for
// preconditioning (here they are the same)
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 5
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(),
             this->matrix.getMat());
#else
  PETSc_call(KSPSetOperators, ksp, this->matrix.getMat(), this->matrix.getMat(),
             SAME_NONZERO_PATTERN);
#endif

  // If this is not called the solution vector is zeroed in the call to
  // KSPSolve().
  PETSc_call(KSPSetInitialGuessNonzero, ksp, PETSC_TRUE);
  PETSc_call(KSPSetFromOptions, ksp);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::solve() {
  Vec & rhs(getDOFManager()._getResidual());
  Vec & solution(getDOFManager()._getSolution());

  this->setOperators();
  MatView(this->matrix.getMat(), PETSC_VIEWER_STDOUT_WORLD);
  VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);

  PETSc_call(KSPSolve, ksp, rhs, solution);
  VecView(solution, PETSC_VIEWER_STDOUT_WORLD);

  this->dof_manager.splitSolutionPerDOFs();
}

/* -------------------------------------------------------------------------- */

DOFManagerPETSc & SparseSolverPETSc::getDOFManager() {
  return aka::as_type<DOFManagerPETSc &>(this->dof_manager);
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::updateInternalParameters() {

  PetscOptionsInsertString(nullptr, petsc_options.c_str());
  KSPSetFromOptions(ksp);
  PetscOptionsClear(nullptr);
}

/* -------------------------------------------------------------------------- */
void SparseSolverPETSc::parseSection(const ParserSection & section) {
  auto parameters = section.getParameters();
  for (auto && param : range(parameters.first, parameters.second)) {
    PetscOptionsSetValue(nullptr, param.getName().c_str(),
                         param.getValue().c_str());
  }
  KSPSetFromOptions(ksp);
  PetscOptionsClear(nullptr);
}

} // namespace akantu
