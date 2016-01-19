/**
 * @file   solver.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Wed Oct 07 2015
 *
 * @brief  Solver interface class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

SolverOptions _solver_no_options(true);

/* -------------------------------------------------------------------------- */
Solver::Solver(SparseMatrix & matrix, const ID & id, const MemoryID & memory_id)
    : Memory(id, memory_id), StaticSolverEventHandler(), matrix(&matrix),
      is_matrix_allocated(false), mesh(NULL),
      communicator(StaticCommunicator::getStaticCommunicator()), solution(NULL),
      synch_registry(NULL) {
  AKANTU_DEBUG_IN();
  StaticSolver::getStaticSolver().registerEventHandler(*this);
  // createSynchronizerRegistry();
  this->synch_registry = new SynchronizerRegistry(*this);
  synch_registry->registerSynchronizer(this->matrix->getDOFSynchronizer(),
                                       _gst_solver_solution);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Solver::~Solver() {
  AKANTU_DEBUG_IN();

  this->destroyInternalData();
  delete synch_registry;
  StaticSolver::getStaticSolver().unregisterEventHandler(*this);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Solver::beforeStaticSolverDestroy() {
  AKANTU_DEBUG_IN();

  try {
    this->destroyInternalData();
  } catch (...) {
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Solver::createSynchronizerRegistry() {
  // this->synch_registry = new SynchronizerRegistry(this);
}

__END_AKANTU__
