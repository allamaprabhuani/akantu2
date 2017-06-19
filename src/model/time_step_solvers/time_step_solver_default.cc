/**
 * @file   time_step_solver_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Sep 16 10:20:55 2015
 *
 * @brief  Default implementation of the time step solver
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
#include "time_step_solver_default.hh"
#include "dof_manager_default.hh"
#include "integration_scheme_1st_order.hh"
#include "integration_scheme_2nd_order.hh"
#include "mesh.hh"
#include "non_linear_solver.hh"
#include "pseudo_time.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::TimeStepSolverDefault(
    DOFManagerDefault & dof_manager, const TimeStepSolverType & type,
    NonLinearSolver & non_linear_solver, const ID & id, UInt memory_id)
    : TimeStepSolver(dof_manager, type, non_linear_solver, id, memory_id),
      dof_manager(dof_manager), is_mass_lumped(false) {
  switch (type) {
  case _tsst_dynamic:
    break;
  case _tsst_dynamic_lumped:
    this->is_mass_lumped = true;
    break;
  case _tsst_static:
    /// initialize a static time solver for callback dofs
    break;
  }
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::setIntegrationScheme(
    const ID & dof_id, const IntegrationSchemeType & type,
    IntegrationScheme::SolutionType solution_type) {
  if (this->integration_schemes.find(dof_id) !=
      this->integration_schemes.end()) {
    AKANTU_EXCEPTION("Their DOFs "
                     << dof_id
                     << "  have already an integration scheme associated");
  }

  IntegrationScheme * integration_scheme = NULL;
  if (this->is_mass_lumped) {
    switch (type) {
    case _ist_forward_euler: {
      integration_scheme = new ForwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_central_difference: {
      integration_scheme = new CentralDifference(dof_manager, dof_id);
      break;
    }
    default:
      AKANTU_EXCEPTION(
          "This integration scheme cannot be used in lumped dynamic");
    }
  } else {
    switch (type) {
    case _ist_pseudo_time: {
      integration_scheme = new PseudoTime(dof_manager, dof_id);
      break;
    }
    case _ist_forward_euler: {
      integration_scheme = new ForwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_trapezoidal_rule_1: {
      integration_scheme = new TrapezoidalRule1(dof_manager, dof_id);
      break;
    }
    case _ist_backward_euler: {
      integration_scheme = new BackwardEuler(dof_manager, dof_id);
      break;
    }
    case _ist_central_difference: {
      integration_scheme = new CentralDifference(dof_manager, dof_id);
      break;
    }
    case _ist_fox_goodwin: {
      integration_scheme = new FoxGoodwin(dof_manager, dof_id);
      break;
    }
    case _ist_trapezoidal_rule_2: {
      integration_scheme = new TrapezoidalRule2(dof_manager, dof_id);
      break;
    }
    case _ist_linear_acceleration: {
      integration_scheme = new LinearAceleration(dof_manager, dof_id);
      break;
    }
    case _ist_generalized_trapezoidal: {
      integration_scheme = new GeneralizedTrapezoidal(dof_manager, dof_id);
      break;
    }
    case _ist_newmark_beta:
      integration_scheme = new NewmarkBeta(dof_manager, dof_id);
      break;
    }
  }

  AKANTU_DEBUG_ASSERT(integration_scheme != nullptr,
                      "No integration scheme was found for the provided types");
  this->integration_schemes[dof_id] = integration_scheme;
  this->solution_types[dof_id] = solution_type;

  this->integration_schemes_owner.insert(dof_id);
}

/* -------------------------------------------------------------------------- */
bool TimeStepSolverDefault::hasIntegrationScheme(const ID & dof_id) const {
  return this->integration_schemes.find(dof_id) !=
         this->integration_schemes.end();
}
/* -------------------------------------------------------------------------- */
TimeStepSolverDefault::~TimeStepSolverDefault() {
  DOFsIntegrationSchemesOwner::iterator it =
      this->integration_schemes_owner.begin();
  DOFsIntegrationSchemesOwner::iterator end =
      this->integration_schemes_owner.end();

  for (; it != end; ++it) {
    delete this->integration_schemes[*it];
  }
  this->integration_schemes.clear();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::solveStep(SolverCallback & solver_callback) {
  this->solver_callback = &solver_callback;

  this->non_linear_solver.solve(*this);

  this->solver_callback = NULL;
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::predictor() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::predictor();

  IntegrationScheme * integration_scheme;
  ID dof_id;

  for(auto & pair : this->integration_schemes) {
    std::tie(dof_id, integration_scheme) = pair;

    if (this->dof_manager.hasPreviousDOFs(dof_id)) {
      this->dof_manager.savePreviousDOFs(dof_id);
    }

    /// integrator predictor
    integration_scheme->predictor(this->time_step);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::corrector() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::corrector();

  IntegrationScheme * integration_scheme;
  ID dof_id;

  for(auto & pair : this->integration_schemes) {
    std::tie(dof_id, integration_scheme) = pair;

    const auto & solution_type = this->solution_types[dof_id];
    integration_scheme->corrector(solution_type, this->time_step);

    /// computing the increment of dof if needed
    if (this->dof_manager.hasDOFsIncrement(dof_id)) {
      if (!this->dof_manager.hasPreviousDOFs(dof_id)) {
        AKANTU_DEBUG_WARNING("In order to compute the increment of "
                             << dof_id << " a 'previous' has to be registered");
        continue;
      }

      Array<Real> & increment = this->dof_manager.getDOFsIncrement(dof_id);
      Array<Real> & previous = this->dof_manager.getPreviousDOFs(dof_id);

      UInt dof_array_comp = this->dof_manager.getDOFs(dof_id).getNbComponent();

      auto prev_dof_it = previous.begin(dof_array_comp);
      auto incr_it = increment.begin(dof_array_comp);
      auto incr_end = increment.end(dof_array_comp);

      increment.copy(this->dof_manager.getDOFs(dof_id));
      for (; incr_it != incr_end; ++incr_it, ++prev_dof_it) {
        *incr_it -= *prev_dof_it;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleJacobian() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleJacobian();

  IntegrationScheme * integration_scheme;
  ID dof_id;

  for(auto & pair : this->integration_schemes) {
    std::tie(dof_id, integration_scheme) = pair;

    const auto & solution_type = this->solution_types[dof_id];

    integration_scheme->assembleJacobian(solution_type,
                                         this->time_step);
  }

  this->dof_manager.applyBoundary("J");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void TimeStepSolverDefault::assembleResidual() {
  AKANTU_DEBUG_IN();

  TimeStepSolver::assembleResidual();

  IntegrationScheme * integration_scheme;
  ID dof_id;

  for(auto & pair : this->integration_schemes) {
    std::tie(dof_id, integration_scheme) = pair;
    integration_scheme->assembleResidual(this->is_mass_lumped);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

} // akantu
