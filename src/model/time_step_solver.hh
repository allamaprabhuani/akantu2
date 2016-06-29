/**
 * @file   time_step_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 24 12:42:04 2015
 *
 * @brief  This corresponding to the time step evolution solver
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
#include "aka_memory.hh"
#include "aka_array.hh"
#include "solver_callback.hh"
#include "integration_scheme.hh"
#include "parameter_registry.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TIME_STEP_SOLVER_HH__
#define __AKANTU_TIME_STEP_SOLVER_HH__

namespace akantu {
class DOFManager;
class NonLinearSolver;
}

__BEGIN_AKANTU__

class TimeStepSolver : public Memory, public ParameterRegistry, public SolverCallback {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  TimeStepSolver(DOFManager & dof_manager, const TimeStepSolverType & type,
                 NonLinearSolver & non_linear_solver, const ID & id,
                 UInt memory_id);
  virtual ~TimeStepSolver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solves on step
  virtual void solveStep(SolverCallback & solver_callback) = 0;

  /// register an integration scheme for a given dof
  virtual void
  setIntegrationScheme(const ID & dof_id, const IntegrationSchemeType & type,
                       IntegrationScheme::SolutionType
                           solution_type = IntegrationScheme::_not_defined) = 0;

  /// replies if a integration scheme has been set
  virtual bool hasIntegrationScheme(const ID & dof_id) const = 0;
  /* ------------------------------------------------------------------------ */
  /* Solver Callback interface                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// implementation of the SolverCallback::predictor()
  virtual void predictor();
  /// implementation of the SolverCallback::corrector()
  virtual void corrector();
  /// implementation of the SolverCallback::assembleJacobian()
  virtual void assembleJacobian();
  /// implementation of the SolverCallback::assembleResidual()
  virtual void assembleResidual();

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(TimeStep, time_step, Real);
  AKANTU_SET_MACRO(TimeStep, time_step, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Underlying dof manager containing the dof to treat
  DOFManager & _dof_manager;

  /// Type of solver
  TimeStepSolverType type;

  /// The time step for this solver
  Real time_step;

  /// Temporary storage for solver callback
  SolverCallback * solver_callback;

  /// NonLinearSolver used by this tome step solver
  NonLinearSolver & non_linear_solver;
};

__END_AKANTU__

#endif /* __AKANTU_TIME_STEP_SOLVER_HH__ */
