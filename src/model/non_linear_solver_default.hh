/**
 * @file   non_linear_solver_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 25 00:48:07 2015
 *
 * @brief Default implementation of NonLinearSolver, in case no external library
 * is there to do the job
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
#include "non_linear_solver.hh"
#include "solver_mumps.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__
#define __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__

namespace akantu {
  class DOFManagerDefault;
}

__BEGIN_AKANTU__

class NonLinearSolverDefault : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverDefault(DOFManagerDefault & dof_manager,
                         const NonLinearSolverType & non_linear_solver_type,
                         const ID & id = "non_linear_solver_default",
                         UInt memory_id = 0);
  virtual ~NonLinearSolverDefault();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Function that solve the non linear system described by the dof manager and
  /// the solver callback functions
  void solve();

protected:
  /// test the convergence compare norm of array to convergence_criteria
  bool testConvergence(const Array<Real> & array);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setParameters(const ParserSection & parameters_section);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  DOFManagerDefault & dof_manager;

  /// Sparse solver used for the linear solves
  SparseSolverMumps solver;

  /// Type of convergence criteria
  SolveConvergenceCriteria convergence_criteria_type;

  /// convergence threshold
  Real convergence_criteria;

  /// Max number of iterations
  UInt max_iterations;

  /// Number of iterations at last solve call
  UInt n_iter;

  /// Convergence error at last solve call
  Real error;

  /// Did the last call to solve reached convergence
  bool converged;
};

__END_AKANTU__

#endif /* __AKANTU_NON_LINEAR_SOLVER_DEFAULT_HH__ */
