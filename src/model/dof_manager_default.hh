/**
 * @file   dof_manager_default.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 11 14:06:18 2015
 *
 * @brief  Default implementation of the dof manager
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
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_DOF_MANAGER_DEFAULT_HH__
#define __AKANTU_DOF_MANAGER_DEFAULT_HH__

__BEGIN_AKANTU__

class SparseMatrixAIJ;

class DOFManagerDefault : public DOFManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DOFManagerDefault(const Mesh & mesh, const ID & id = "dof_manager_default",
                    const MemoryID & memory_id = 0);
  virtual ~DOFManagerDefault();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Get the part of the solution corresponding to the dof_id
  virtual void getSolution(const ID & dof_id, Array<Real> & solution_array);

  /// Assemble an array to the global residual array
  virtual void assembleToResidual(const ID & dof_id,
                                  const Array<Real> & array_to_assemble,
                                  Real scale_factor = 1.);

  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalArrayResidual(
      const ID & dof_id, const Array<Real> & array_to_assemble,
      const ElementType & type, const GhostType & ghost_type,
      Real scale_factor = 1.);
  /**
   * Assemble elementary values to the global residual array. The dof number is
   * implicitly considered as conn(el, n) * nb_nodes_per_element + d.
   * With 0 < n < nb_nodes_per_element and 0 < d < nb_dof_per_node
   **/
  virtual void assembleElementalMatricesToMatrix(
      const ID & matrix_id, const ID & dof_id,
      const Array<Real> & elementary_mat, const ElementType & type,
      const GhostType & ghost_type, const MatrixType & elemental_matrix_type,
      const Array<UInt> & filter_elements);

private:
  /// Add a symmetric matrices to a symmetric sparse matrix
  void addSymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a unsymmetric matrices to a symmetric sparse matrix (i.e. cohesive
  /// elements)
  void addUnsymmetricElementalMatrixToSymmetric(
      SparseMatrixAIJ & matrix, const Matrix<Real> & element_mat,
      const Vector<UInt> & equation_numbers, UInt max_size);

  /// Add a matrices to a unsymmetric sparse matrix
  void addElementalMatrixToUnsymmetric(SparseMatrixAIJ & matrix,
                                       const Matrix<Real> & element_mat,
                                       const Vector<UInt> & equation_numbers,
                                       UInt max_size);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get an instance of a new SparseMatrix
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const MatrixType & matrix_type);

  /// Get an instance of a new SparseMatrix as a copy of the SparseMatrix
  /// matrix_to_copy_id
  virtual SparseMatrix & getNewMatrix(const ID & matrix_id,
                                      const ID & matrix_to_copy_id);

  /// Get the reference of an existing matrix
  SparseMatrixAIJ & getMatrix(const ID & matrix_id);

  /// Get the solution array
  AKANTU_GET_MACRO_NOT_CONST(Solution, solution, Array<Real> &);
  /// Get the residual array
  AKANTU_GET_MACRO_NOT_CONST(Residual, residual, Array<Real> &);
  /// Get the blocked dofs array
  AKANTU_GET_MACRO(BlockedDOFs, blocked_dofs, const Array<bool> &);

  bool isLocalOrMasterDOF(UInt dof_num);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<ID, SparseMatrixAIJ *> AIJMatrixMap;

  /// rhs to the system of equation corresponding to the residual linked to the
  /// different dofs
  Array<Real> residual;

  /// solution of the system of equation corresponding to the different dofs
  Array<Real> solution;

  /// blocked degree of freedom in the system equation corresponding to the
  /// different dofs
  Array<bool> blocked_dofs;

  /// define the dofs type, local, shared, ghost
  Array<Int> dofs_type;

  /// Map of the different matrices stored in the dof manager
  AIJMatrixMap aij_matrices;
};

__END_AKANTU__

#include "dof_manager_default_inline_impl.cc"

#endif /* __AKANTU_DOF_MANAGER_DEFAULT_HH__ */
