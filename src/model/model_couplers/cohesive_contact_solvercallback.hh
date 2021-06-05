
/**
 * @file   cohesive_contact_solvercallback.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu Jan 17 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model via solvercallback
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
#include "mesh_iterators.hh"
#include "non_linear_solver.hh"
#include "contact_mechanics_model.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */


#ifndef __AKANTU_COHESIVE_CONTACT_SOLVERCALLBACK_HH__
#define __AKANTU_COHESIVE_CONTACT_SOLVERCALLBACK_HH__

namespace akantu{


class CohesiveContactSolverCallback : public SolverCallback {

public:
  CohesiveContactSolverCallback(SolidMechanicsModelCohesive &,
                                ContactMechanicsModel &, AnalysisMethod &);

public:
  void assembleMatrix(const ID &) override;

  void assembleResidual() override;

  void assembleLumpedMatrix(const ID &) override;

  MatrixType getMatrixType(const ID &) override;

  void predictor() override;

  void corrector() override;

  void beforeSolveStep() override;

  void afterSolveStep(bool converged=true) override;

private:
  SolidMechanicsModelCohesive &solid;

  ContactMechanicsModel &contact;

  AnalysisMethod & method;
};

}
  
#endif
