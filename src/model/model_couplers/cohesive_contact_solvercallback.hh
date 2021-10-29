/**
 * @file   cohesive_contact_solvercallback.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Jun 06 2021
 * @date last modification: Wed Jul 28 2021
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model via solvercallback
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "contact_mechanics_model.hh"
#include "mesh_iterators.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_CONTACT_SOLVERCALLBACK_HH__
#define __AKANTU_COHESIVE_CONTACT_SOLVERCALLBACK_HH__

namespace akantu {

class CohesiveContactSolverCallback : public SolverCallback {

public:
  CohesiveContactSolverCallback(SolidMechanicsModelCohesive & /*solid*/,
                                ContactMechanicsModel & /*contact*/,
                                AnalysisMethod & /*method*/);

public:
  /// implementation of SolverCallback::assembleMatrix
  void assembleMatrix(const ID & /*matrix_id*/) override;

  /// implementation of SolverCallback::assembleResidual
  void assembleResidual() override;

  /// implementation of SolverCallback::assembleLumpedMatrix
  void assembleLumpedMatrix(const ID & /*matrix_id*/) override;

  /// implementation of SolverCallback::getMatrixType
  MatrixType getMatrixType(const ID & /*unused*/) override;

  /// implementation of SolverCallback::predictor
  void predictor() override;

  /// implementation of SolverCallback::corrector
  void corrector() override;

  /// implementation of SolverCallback::beforeSolveStep
  void beforeSolveStep() override;

  /// implementation of SolverCallback::afterSolveStep
  void afterSolveStep(bool converged = true) override;

private:
  /// model for the solid mechanics part of the coupling
  SolidMechanicsModelCohesive & solid;

  /// model for the contact resoluion of the coupling
  ContactMechanicsModel & contact;

  /// Method of resolution for the coupling solver
  AnalysisMethod & method;
};

} // namespace akantu

#endif
