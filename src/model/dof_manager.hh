/**
 * @file   dof_manager.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Jul 22 11:43:43 2015
 *
 * @brief  Class handling the different types of dofs
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
#include "non_linear_solver.hh"
#include "time_step_solver.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_HH__
#define __AKANTU_DOF_MANAGER_HH__

namespace akantu {
class SolverCallback;
}

__BEGIN_AKANTU__

class DOFManager : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManager(const Mesh & mesh, SolverCallback & solver_callback,
             const ID & id = "dof_manager", const MemoryID & memory_id = 0);
  virtual ~DOFManager();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// register an array of degree of freedom
  void registerDOFs(const ID & dof_id, Array<Real> & dofs_array);

  /// register an array of derivatives for a particular dof array
  void registerDOFDerivative(const ID & dof_id, UInt order,
                             Array<Real> & dofs_derivative);

  /// register array representing the blocked degree of freedoms
  void registerBlockedDOFs(const ID & dof_id, Array<Real> & blocked_dofs);

  /// Get the part of the solution corresponding to the dof_id
  virtual void getSolution(const ID & dof_id, Array<Real> & solution_array) = 0;

  /// Assemble an array to the global residual array
  virtual void assembleToResidual(const ID & dof_id,
                                  const Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.) = 0;

  /**
   * Assemble elementary values to a local array of the size nb_nodes *
   * nb_dof_per_node. The dof number is implicitly considered as
   * conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayLocalArray(
      const Array<Real> & elementary_vect, Array<Real> & array_assembeled,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayResidual(
      const ID & dof_id, const Array<Real> & elementary_vect,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.,
      const Array<UInt> & filter_elements = empty_filter);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.  With 0 <
   * n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, const ElementType & type,
      const GhostType & ghost_type,
      const MatrixType & elemental_matrix_type = _symmetric,
      const Array<UInt> & filter_elements = empty_filter) = 0;

  /// notation fully defined yet...
  // virtual void assemblePreassembledMatrix(const ID & matrix_id,
  //                                         const ID & dof_id_m,
  //                                         const ID & dof_id_n,
  //                                         const Matrix<Real> & matrix) = 0;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// fill a Vector with the equation numbers corresponding to the given
  /// connectivity
  inline void extractElementEquationNumber(
      const Array<UInt> & equation_numbers, const Vector<UInt> & connectivity,
      UInt nb_degree_of_freedom, Vector<UInt> & local_equation_number);

  /// register a matrix
  void registerSparseMatrix(const ID & matrix_id, SparseMatrix & matrix);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the equation numbers (in local numbering) corresponding to a dof ID
  Array<UInt> & getLocalEquationNumbers(const ID & dof_id);
  /// get the equation numbers (in global numbering) corresponding to a dof ID
  Array<UInt> & getGlobalEquationNumbers(const ID & dof_id);

  /// get the equation numbers (in local numbering) corresponding to a dof ID
  const Array<UInt> & getLocalEquationNumbers(const ID & dof_id) const;
  /// get the equation numbers (in global numbering) corresponding to a dof ID
  const Array<UInt> & getGlobalEquationNumbers(const ID & dof_id) const;

  /// Global number of dofs
  AKANTU_GET_MACRO(SystemSize, this->system_size, UInt);

  /// Get the solver callback stored in the dof_manager
  AKANTU_GET_MACRO(SolverCallback, solver_callback, SolverCallback &);

  /* ------------------------------------------------------------------------ */
  /* Matrices accessors                                                       */
  /* ------------------------------------------------------------------------ */
  /// Get an instance of a new SparseMatrix
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const MatrixType & matrix_type);

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const ID & matrix_to_copy_id);

  /// Get the reference of an existing matrix
  SparseMatrix & getMatrix(const ID & matrix_id);

  /* ------------------------------------------------------------------------ */
  /* Non linear system solver                                                 */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a non linear solver
  virtual NonLinearSolver &
  getNewNonLinearSolver(const ID & nls_solver_id,
                        const NonLinearSolverType & _non_linear_solver_type);

  /// get instance of a non linear solver
  virtual NonLinearSolver & getNonLinearSolver(const ID & nls_solver_id);

  /* ------------------------------------------------------------------------ */
  /* TimeStepSolver                                                           */
  /* ------------------------------------------------------------------------ */
  /// Get instance of a time step solver
  virtual TimeStepSolver &
  getNewTimeStepSolver(const ID & time_step_solver_id, UInt order,
                       const TimeStepSolverType & type);

  /// get instance of a time step solver
  virtual TimeStepSolver & getTimeStepSolver(const ID & time_step_solver_id);

  /* ------------------------------------------------------------------------*/
  /* Class Members                                                           */
  /* ------------------------------------------------------------------------*/
protected:
  /// dof representations in the dof manager
  struct DOFData {
    /// Degree of freedom array
    Array<Real> * dof;
    /// Blocked degree of freedoms array
    Array<Real> * blocked_dofs;
    /* ---------------------------------------------------------------------- */
    /* data for dynamic simulations                                           */
    /* ---------------------------------------------------------------------- */
    /// Degree of freedom derivatives arrays
    std::vector<Array<Real> *> dof_derivatives;

    /// global numbering equation numbers
    Array<UInt> global_equation_number;
    /// local numbering equation numbers
    Array<UInt> local_equation_number;
    /// link between global and local equation numbers
    unordered_map<UInt, UInt>::type global_to_local_equation_number;
  };

  typedef std::map<ID, DOFData *> DOFStorage;

  /// store a reference to the dof arrays
  DOFStorage dofs;

  /// list of sparse matrices that where created
  std::map<ID, SparseMatrix *> matrices;

  /// non linear solvers storage
  std::map<ID, NonLinearSolver *> non_linear_solvers;

  /// time step solvers storage
  std::map<ID, TimeStepSolver *> time_step_solvers;

  /// reference to the underlying mesh
  const Mesh & mesh;

  /// Solver callback to assemble residual and jacobian
  SolverCallback & solver_callback;

  /// Total number of degrees of freedom
  UInt local_system_size;

  /// Total number of degrees of freedom
  UInt system_size;
};

__END_AKANTU__

#include "dof_manager_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_HH__ */
