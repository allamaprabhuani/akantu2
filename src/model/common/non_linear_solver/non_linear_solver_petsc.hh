/**
 * @file   non_linear_solver_petsc.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  Tue Jan 01 2019
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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include <petscsnes.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__
#define __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__

namespace akantu {
class DOFManagerPETSc;
} // namespace akanatu

namespace akantu {

class NonLinearSolverPETSc : public NonLinearSolver {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLinearSolverPETSc(DOFManagerPETSc & dof_manager,
                       const NonLinearSolverType & non_linear_solver_type,
                       const ID & id = "non_linear_solver_petsc",
                       UInt memory_id = 0);

  ~NonLinearSolverPETSc() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// solve the system described by the jacobian matrix, and rhs contained in
  /// the dof manager
  void solve(SolverCallback & callback) override;

  /// parse the arguments from the input file
  void parseSection(const ParserSection & section) override;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  DOFManagerPETSc & dof_manager;

  SNES snes;
};

} // namepsace akantu

#endif /* __AKANTU_NON_LINEAR_SOLVER_PETSC_HH__ */
