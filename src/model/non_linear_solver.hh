/**
 * @file   non_linear_solver.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 24 23:48:41 2015
 *
 * @brief  Non linear solver interface
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
#include "aka_memory.hh"
#include "parameter_registry.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_HH__
#define __AKANTU_NON_LINEAR_SOLVER_HH__

namespace akantu {
class DOFManager;
class SolverCallback;
}

__BEGIN_AKANTU__

class NonLinearSolver : Memory, public ParameterRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolver(DOFManager & dof_manager,
                  const NonLinearSolverType & non_linear_solver_type,
                  const ID & id = "non_linear_solver", UInt memory_id = 0);
  virtual ~NonLinearSolver();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  virtual void solve(SolverCallback & callback) = 0;

protected:
  void checkIfTypeIsSupported();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  DOFManager & _dof_manager;

protected:
  /// type of non linear solver
  NonLinearSolverType non_linear_solver_type;

  /// list of supported non linear solver types
  std::set<NonLinearSolverType> supported_type;
};

namespace debug {
  class NLSNotConvergedException : public Exception {
  public:
    NLSNotConvergedException(Real threshold, UInt niter)
        : Exception("The non linear solver did not converge."),
          threshold(threshold), niter(niter) {}
    Real threshold;
    UInt niter;
  };
}

__END_AKANTU__

#endif /* __AKANTU_NON_LINEAR_SOLVER_HH__ */
