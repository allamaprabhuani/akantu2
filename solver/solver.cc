/**
 * @file   solver.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 16:19:27 2010
 *
 * @brief  Solver interface class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solver.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Solver::Solver(SparseMatrix & matrix,
	       const SolverID & id,
	       const MemoryID & memory_id) :
  Memory(memory_id), id(id), matrix(&matrix), is_matrix_allocated(false), mesh(NULL) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
Solver::Solver(const Mesh & mesh,
	       __attribute__ ((unused)) const SparseMatrixType & sparse_matrix_type,
	       __attribute__ ((unused)) UInt nb_degre_of_freedom,
	       const SolverID & id,
	       const MemoryID & memory_id) :
  Memory(memory_id), id(id), matrix(NULL), is_matrix_allocated(false), rhs(NULL), mesh(&mesh) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Solver::~Solver() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
