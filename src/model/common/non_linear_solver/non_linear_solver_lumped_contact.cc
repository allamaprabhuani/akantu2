/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_linear_solver_lumped_contact.hh"
#include "dof_manager_default.hh"
#include "solver_callback.hh"
#include "solver_vector_default.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLinearSolverLumpedContact::NonLinearSolverLumpedContact(
    DOFManagerDefault & dof_manager,
    const NonLinearSolverID & non_linear_solver_type, const ID & id)
    : NonLinearSolverLumped(dof_manager, non_linear_solver_type, id) {
  this->supported_type.insert("lumped_contact");
  this->checkIfTypeIsSupported();
}

/* ------------------------------------------------------------------------ */
void NonLinearSolverLumpedContact::solve(SolverCallback & solver_callback) {
  // --------------------
  // computes everything without contact forces
  solver_callback.beforeSolveStep();
  this->dof_manager.updateGlobalBlockedDofs();
  solver_callback.predictor();

  solver_callback.assembleResidual();

  auto & x = aka::as_type<SolverVectorDefault>(this->dof_manager.getSolution());
  const auto & b = this->dof_manager.getResidual();

  x.resize();

  const auto & blocked_dofs = this->dof_manager.getBlockedDOFs();
  const auto & A = this->dof_manager.getLumpedMatrix("M");

  // alpha is the conversion factor from from force/mass to acceleration needed
  // in model coupled with atomistic \todo find a way to define alpha per dof
  // type
  x.zero();
  NonLinearSolverLumped::solveLumped(A, x, b, alpha, blocked_dofs);

  this->dof_manager.splitSolutionPerDOFs();

  solver_callback.corrector();
  // --------------------

  // computes contact forces
  solver_callback.assembleResidual("contact");

  // restores everything
  solver_callback.afterSolveStep(false);

  // --------------------
  // redoes the step with the contact forces added
  solver_callback.beforeSolveStep();

  solver_callback.predictor();

  x.zero();
  NonLinearSolverLumped::solveLumped(A, x, b, alpha, blocked_dofs);

  this->dof_manager.splitSolutionPerDOFs();
  solver_callback.corrector();

  solver_callback.afterSolveStep(true);
}

/* -------------------------------------------------------------------------- */
[[maybe_unused]] bool non_linear_solver_is_allocated_lumped_contact =
    instantiateNonLinearSolver<NonLinearSolverLumpedContact, DOFManagerDefault>(
        "lumped_contact");

} // namespace akantu
