/**
 * @file   solver.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 16:19:27 2010
 *
 * @brief  Solver interface class
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
#include "solver.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

SolverOptions _solver_no_options(true);

/* -------------------------------------------------------------------------- */
Solver::Solver(SparseMatrix & matrix,
	       const ID & id,
	       const MemoryID & memory_id) :
  Memory(memory_id), id(id), matrix(&matrix), is_matrix_allocated(false), mesh(NULL) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Solver::~Solver() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
