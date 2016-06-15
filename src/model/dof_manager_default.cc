/**
 * @file   dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Aug 11 16:21:01 2015
 *
 * @brief  Implementation of the default DOFManager
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
#include "dof_manager_default.hh"
#include "sparse_matrix_aij.hh"
#include "time_step_solver_default.hh"
#include "static_communicator.hh"
#include "non_linear_solver_default.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addSymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = i; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addUnsymmetricElementalMatrixToSymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          if (c_jcn >= c_irn) {
            matrix(c_irn, c_jcn) += elementary_mat(i, j);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void DOFManagerDefault::addElementalMatrixToUnsymmetric(
    SparseMatrixAIJ & matrix, const Matrix<Real> & elementary_mat,
    const Vector<UInt> & equation_numbers, UInt max_size) {
  for (UInt i = 0; i < elementary_mat.rows(); ++i) {
    UInt c_irn = equation_numbers(i);
    if (c_irn < max_size) {
      for (UInt j = 0; j < elementary_mat.cols(); ++j) {
        UInt c_jcn = equation_numbers(j);
        if (c_jcn < max_size) {
          matrix(c_irn, c_jcn) += elementary_mat(i, j);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::DOFManagerDefault(const ID & id, const MemoryID & memory_id)
    : DOFManager(id, memory_id), residual(0, 1, std::string(id + ":residual")),
      global_solution(0, 1, std::string(id + ":global_solution")),
      global_blocked_dofs(0, 1, std::string(id + ":global_blocked_dofs")),
      dofs_type(0, 1, std::string(id + ":dofs_type")),
      data_cache(0, 1, std::string(id + ":data_cache_array"){}

/* -------------------------------------------------------------------------- */
DOFManagerDefault::~DOFManagerDefault() {}

/* -------------------------------------------------------------------------- */
template <typename T>
void DOFManagerDefault::assembleToGlobalArray(
    const ID & dof_id, const Array<T> & array_to_assemble,
    Array<T> & global_array, T scale_factor) {
  AKANTU_DEBUG_IN();
  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms =
      array_to_assemble.getSize() * array_to_assemble.getNbComponent();

  AKANTU_DEBUG_ASSERT(equation_number.getSize() == nb_degree_of_freedoms,
                      "The array to assemble does not have a correct size."
                          << " (" << array_to_assemble.getID() << ")");

  typename Array<T>::const_scalar_iterator arr_it =
      array_to_assemble.begin_reinterpret(nb_degree_of_freedoms);
  Array<UInt>::const_scalar_iterator equ_it = equation_number.begin();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++arr_it, ++equ_it) {
    global_array(*equ_it) += scale_factor * (*arr_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::registerDOFs(const ID & dof_id,
                                     Array<Real> & dofs_array,
                                     const DOFSupportType & support_type) {
  // stores the current numbers of dofs
  UInt local_nb_dofs = this->local_system_size;
  UInt pure_local_nb_dofs = this->pure_local_system_size;

  // update or create the dof_data
  DOFManager::registerDOFs(dof_id, dofs_array, support_type);

  // Count the number of pure local dofs per proc
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  UInt prank = comm.whoAmI();
  UInt psize = comm.getNbProc();

  Array<UInt> nb_dofs_per_proc(psize);
  nb_dofs_per_proc(prank) = this->pure_local_system_size - pure_local_nb_dofs;
  comm.allGather(nb_dofs_per_proc);

  UInt first_global_dofs_id = std::accumulate(
      nb_dofs_per_proc.begin(), nb_dofs_per_proc.begin() + prank, 0);

  // nb local dofs to account for
  UInt nb_dofs = this->local_system_size - local_nb_dofs;

  DOFData & dof_data = *dofs[dof_id];

  this->global_equation_number.resize(this->local_system_size);

  // set the equation numbers
  UInt first_dof_id = local_nb_dofs;
  dof_data.local_equation_number.resize(nb_dofs);
  this->dofs_type.resize(local_system_size);

  for (UInt d = 0; d < nb_dofs; ++d) {
    UInt local_eq_num = first_dof_id + d;
    dof_data.local_equation_number(d) = local_eq_num;
    UInt global_eq_num = first_global_dofs_id + d;
    this->global_equation_number(local_eq_num) = global_eq_num;
    this->global_to_local_mapping[global_eq_num] = local_eq_num;
    switch (support_type) {
    case _dst_nodal: {
      UInt node = d / dof_data.dof->getNbComponent();
      this->dofs_type(local_eq_num) = this->mesh->getNodeType(node);
      break;
    }
    case _dst_generic: {
      this->dofs_type(local_eq_num) = _nt_normal;
      break;
    }
    default: { AKANTU_EXCEPTION("This type of dofs is not handled yet."); }
    }
  }

  this->residual.resize(this->local_system_size);
  this->global_solution.resize(this->local_system_size);
  this->global_blocked_dofs.resize(this->local_system_size);
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const MatrixType & matrix_type) {
  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatrix * sm = new SparseMatrixAIJ(*this, matrix_type, matrix_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrix & DOFManagerDefault::getNewMatrix(const ID & id,
                                               const ID & matrix_to_copy_id) {

  ID matrix_id = this->id + ":mtx:" + id;
  SparseMatrixAIJ & sm_to_copy = this->getMatrix(matrix_to_copy_id);
  SparseMatrix * sm = new SparseMatrixAIJ(sm_to_copy, matrix_id);
  this->registerSparseMatrix(matrix_id, *sm);

  return *sm;
}

/* -------------------------------------------------------------------------- */
SparseMatrixAIJ & DOFManagerDefault::getMatrix(const ID & id) {
  SparseMatrix & matrix = DOFManager::getMatrix(id);

  return dynamic_cast<SparseMatrixAIJ &>(matrix);
}

/* -------------------------------------------------------------------------- */
NonLinearSolver &
DOFManagerDefault::getNewNonLinearSolver(const ID & id,
                                         const NonLinearSolverType & type) {
  ID non_linear_solver_id = this->id + ":nls:" + id;
  NonLinearSolver * nls = NULL;
  switch (type) {
  case _nls_newton_raphson:
  case _nls_newton_raphson_modified: {
    nls = new NonLinearSolverNewtonRaphson(*this, type, non_linear_solver_id,
                                           this->memory_id);
    break;
  }
  case _nls_linear: {
    nls = new NonLinearSolverLinear(*this, type, non_linear_solver_id,
                                    this->memory_id);
    break;
  }
  case _nls_lumped: {
    nls = new NonLinearSolverLumped(*this, type, non_linear_solver_id,
                                    this->memory_id);
    break;
  }
  default:
    AKANTU_EXCEPTION("The asked type of non linear solver is not supported by "
                     "this dof manager");
  }

  this->registerNonLinearSolver(non_linear_solver_id, *nls);

  return *nls;
}

/* -------------------------------------------------------------------------- */
TimeStepSolver &
DOFManagerDefault::getNewTimeStepSolver(const ID & id,
                                        const TimeStepSolverType & type,
                                        NonLinearSolver & non_linear_solver) {
  ID time_step_solver_id = this->id + ":tss:" + id;

  TimeStepSolver * tss = new TimeStepSolverDefault(
      *this, type, non_linear_solver, time_step_solver_id, this->memory_id);

  this->registerTimeStepSolver(time_step_solver_id, *tss);

  return *tss;
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearResidual() {
  this->residual.resize(this->local_system_size);
  this->residual.clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearMatrix(const ID & mtx) {
  this->getMatrix(mtx).clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::clearLumpedMatrix(const ID & mtx) {
  this->getLumpedMatrix(mtx).clear();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::updateGlobalBlockedDofs() {
  DOFStorage::iterator it = this->dofs.begin();
  DOFStorage::iterator end = this->dofs.end();

  this->global_blocked_dofs.resize(this->local_system_size);
  this->global_blocked_dofs.clear();

  for (; it != end; ++it) {
    DOFData & dof_data = *it->second;
    this->assembleToGlobalArray(it->first, *dof_data.blocked_dofs,
                                this->global_blocked_dofs, true);
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getArrayPerDOFs(const ID & dof_id,
                                        const Array<Real> & global_array,
                                        Array<Real> & local_array) const {
  AKANTU_DEBUG_IN();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  UInt nb_degree_of_freedoms = equation_number.getSize();
  local_array.resize(nb_degree_of_freedoms / local_array.getNbComponent());

  Array<Real>::scalar_iterator loc_it =
      local_array.begin_reinterpret(nb_degree_of_freedoms);
  Array<UInt>::const_scalar_iterator equ_it = equation_number.begin();

  for (UInt d = 0; d < nb_degree_of_freedoms; ++d, ++loc_it, ++equ_it) {
    (*loc_it) = global_array(*equ_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getSolutionPerDOFs(const ID & dof_id,
                                           Array<Real> & solution_array) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->global_solution, solution_array);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
void DOFManagerDefault::getLumpedMatrixPerDOFs(const ID & dof_id,
                                               const ID & lumped_mtx,
                                               Array<Real> & lumped) {
  AKANTU_DEBUG_IN();
  this->getArrayPerDOFs(dof_id, this->getLumpedMatrix(lumped_mtx), lumped);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleToResidual(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    Real scale_factor) {
  AKANTU_DEBUG_IN();

  this->assembleToGlobalArray(dof_id, array_to_assemble, this->residual,
                              scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleToLumpedMatrix(
    const ID & dof_id, const Array<Real> & array_to_assemble,
    const ID & lumped_mtx, Real scale_factor) {
  AKANTU_DEBUG_IN();

  Array<Real> & lumped = this->getLumpedMatrix(lumped_mtx);
  this->assembleToGlobalArray(dof_id, array_to_assemble, lumped, scale_factor);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleMatMulVectToResidual(const ID & dof_id,
                                                     const ID & A_id,
                                                     const Array<Real> & x,
                                                     Real scale_factor) {
  SparseMatrixAIJ & A = this->getMatrix(A_id);
  this->data_cache.resize(this->local_system_size);
  this->data_cache.clear();

  this->assembleToGlobalArray(dof_id, x, this->data_cache, 1.);

  A.matVecMul(this->data_cache, this->residual, scale_factor, 1.);
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleLumpedMatMulVectToResidual(const ID & dof_id,
                                                           const ID & A_id,
                                                           const Array<Real> & x,
                                                           Real scale_factor) {
  const Array<Real> & A = this->getLumpedMatrix(A_id);

  this->data_cache.resize(this->local_system_size);
  this->data_cache.clear();
  this->assembleToGlobalArray(dof_id, x, this->data_cache, scale_factor);

  Array<Real>::const_scalar_iterator A_it = A.begin();
  Array<Real>::const_scalar_iterator A_end = A.end();
  Array<Real>::const_scalar_iterator x_it = this->data_cache.begin();
  Array<Real>::scalar_iterator r_it = this->residual.begin();

  for (; A_it != A_end; ++A_it, ++x_it, ++r_it) {
    *r_it += *A_it ** x_it;
  }
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::assembleElementalMatricesToMatrix(
    const ID & matrix_id, const ID & dof_id, const Array<Real> & elementary_mat,
    const ElementType & type, const GhostType & ghost_type,
    const MatrixType & elemental_matrix_type,
    const Array<UInt> & filter_elements) {
  AKANTU_DEBUG_IN();

  this->addToProfile(matrix_id, dof_id);

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);
  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  UInt nb_element;
  if (ghost_type == _not_ghost) {
    nb_element = this->mesh->getNbElement(type);
  } else {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

  UInt * filter_it = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filter_it = filter_elements.storage();
  } else {
    nb_element = this->mesh->getNbElement(type, ghost_type);
  }

  AKANTU_DEBUG_ASSERT(elementary_mat.getSize() == nb_element,
                      "The vector elementary_mat("
                          << elementary_mat.getID()
                          << ") has not the good size.");

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

  UInt nb_degree_of_freedom = this->getDOFs(dof_id).getNbComponent();
  // UInt nb_degree_of_freedom = elementary_mat.getNbComponent() /
  //                             (nb_nodes_per_element * nb_nodes_per_element);

  const Array<UInt> & connectivity =
      this->mesh->getConnectivity(type, ghost_type);
  Array<UInt>::const_vector_iterator conn_begin =
      connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator conn_it = conn_begin;

  UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom;

  Vector<UInt> element_eq_nb(nb_degree_of_freedom * nb_nodes_per_element);
  Array<Real>::const_matrix_iterator el_mat_it =
      elementary_mat.begin(size_mat, size_mat);

  for (UInt e = 0; e < nb_element; ++e, ++el_mat_it) {
    if (filter_it != NULL)
      conn_it = conn_begin + *filter_it;

    this->extractElementEquationNumber(equation_number, *conn_it,
                                       nb_degree_of_freedom, element_eq_nb);
    this->localToGlobalEquationNumber(element_eq_nb);

    if (filter_it != NULL)
      ++filter_it;
    else
      ++conn_it;

    if (A.getMatrixType() == _symmetric)
      if (elemental_matrix_type == _symmetric)
        this->addSymmetricElementalMatrixToSymmetric(
            A, *el_mat_it, element_eq_nb, A.getSize());
      else
        this->addUnsymmetricElementalMatrixToSymmetric(
            A, *el_mat_it, element_eq_nb, A.getSize());
    else
      this->addElementalMatrixToUnsymmetric(A, *el_mat_it, element_eq_nb,
                                            A.getSize());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void DOFManagerDefault::addToProfile(const ID & matrix_id, const ID & dof_id) {
  AKANTU_DEBUG_IN();

  const DOFData & dof_data = this->getDOFData(dof_id);

  if (dof_data.support_type != _dst_nodal)
    return;

  std::pair<ID, ID> mat_dof = std::make_pair(matrix_id, dof_id);
  if (this->matrix_profiled_dofs.find(mat_dof) !=
      this->matrix_profiled_dofs.end())
    return;

  UInt nb_degree_of_freedom_per_node = dof_data.dof->getNbComponent();

  const Array<UInt> & equation_number = this->getLocalEquationNumbers(dof_id);

  SparseMatrixAIJ & A = this->getMatrix(matrix_id);

  // if(irn_jcn_to_k) delete irn_jcn_to_k;
  // irn_jcn_to_k = new std::map<std::pair<UInt, UInt>, UInt>;
  //  A.clearProfile();

  UInt size = A.getSize();

  Mesh::type_iterator it = this->mesh->firstType(
      this->mesh->getSpatialDimension(), _not_ghost, _ek_not_defined);
  Mesh::type_iterator end = this->mesh->lastType(
      this->mesh->getSpatialDimension(), _not_ghost, _ek_not_defined);
  for (; it != end; ++it) {
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

    const Array<UInt> & connectivity =
        this->mesh->getConnectivity(*it, _not_ghost);
    Array<UInt>::const_vector_iterator cit =
        connectivity.begin(nb_nodes_per_element);
    Array<UInt>::const_vector_iterator cend =
        connectivity.end(nb_nodes_per_element);

    UInt size_mat = nb_nodes_per_element * nb_degree_of_freedom_per_node;
    Vector<UInt> element_eq_nb(size_mat);

    for (; cit != cend; ++cit) {
      this->extractElementEquationNumber(
          equation_number, *cit, nb_degree_of_freedom_per_node, element_eq_nb);
      this->localToGlobalEquationNumber(element_eq_nb);

      for (UInt i = 0; i < size_mat; ++i) {
        UInt c_irn = element_eq_nb(i);
        if (c_irn < size) {
          for (UInt j = 0; j < size_mat; ++j) {
            UInt c_jcn = element_eq_nb(j);
            if (c_jcn < size) {
              A.addToProfile(c_irn, c_jcn);
            }
          }
        }
      }
    }
  }

  this->matrix_profiled_dofs.insert(mat_dof);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
