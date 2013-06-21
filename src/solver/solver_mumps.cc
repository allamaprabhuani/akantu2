/**
 * @file   solver_mumps.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Dec 13 10:48:06 2010
 *
 * @brief  implem of SolverMumps class
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
 * @section DESCRIPTION
 *
 * @subsection Ctrl_param Control parameters
 *
 * ICNTL(1),
 * ICNTL(2),
 * ICNTL(3) : output streams for error, diagnostics, and global messages
 *
 * ICNTL(4) : verbose level : 0 no message - 4 all messages
 *
 * ICNTL(5) : type of matrix, 0 assembled, 1 elementary
 *
 * ICNTL(6) : control  the permutation and scaling(default 7)  see mumps doc for
 * more information
 *
 * ICNTL(7) : determine  the pivot  order (default  7) see  mumps doc  for more
 * information
 *
 * ICNTL(8) : describe the scaling method used
 *
 * ICNTL(9) : 1 solve A x = b, 0 solve At x = b
 *
 * ICNTL(10) : number of iterative refinement when NRHS = 1
 *
 * ICNTL(11) : > 0 return statistics
 *
 * ICNTL(12) : only used for SYM = 2, ordering strategy
 *
 * ICNTL(13) :
 *
 * ICNTL(14) : percentage of increase of the estimated working space
 *
 * ICNTL(15-17) : not used
 *
 * ICNTL(18) : only  used if ICNTL(5) = 0, 0 matrix  centralized, 1 structure on
 * host and mumps  give the mapping, 2 structure on  host and distributed matrix
 * for facto, 3 distributed matrix
 *
 * ICNTL(19) : > 0, Shur complement returned
 *
 * ICNTL(20) : 0 rhs dense, 1 rhs sparse
 *
 * ICNTL(21) : 0 solution in rhs, 1 solution distributed in ISOL_loc and SOL_loc
 * allocated by user
 *
 * ICNTL(22) : 0 in-core, 1 out-of-core
 *
 * ICNTL(23) : maximum memory allocatable by mumps pre proc
 *
 * ICNTL(24) : controls the detection of "null pivot rows"
 *
 * ICNTL(25) :
 *
 * ICNTL(26) :
 *
 * ICNTL(27) :
 *
 * ICNTL(28) : 0 automatic choice, 1 sequential analysis, 2 parallel analysis
 *
 * ICNTL(29) : 0 automatic choice, 1 PT-Scotch, 2 ParMetis
 */

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_MPI
#include "static_communicator_mpi.hh"
#endif

#include "solver_mumps.hh"
#include "dof_synchronizer.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolverMumps::SolverMumps(SparseMatrix & matrix,
			 const ID & id,
			 const MemoryID & memory_id) :
  Solver(matrix, id, memory_id), is_mumps_data_initialized(false), rhs_is_local(true) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_MPI
  parallel_method = SolverMumpsOptions::_fully_distributed;
#else //AKANTU_USE_MPI
  parallel_method = SolverMumpsOptions::_master_slave_distributed;
#endif //AKANTU_USE_MPI

  CommunicatorEventHandler & comm_event_handler = *this;

  communicator.registerEventHandler(comm_event_handler);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::~SolverMumps() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::destroyMumpsData() {
  AKANTU_DEBUG_IN();

  if(is_mumps_data_initialized) {
    mumps_data.job = _smj_destroy; // destroy
    dmumps_c(&mumps_data);
    is_mumps_data_initialized = false;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::onCommunicatorFinalize(const StaticCommunicator & comm) {
  AKANTU_DEBUG_IN();

  try{
#if defined(AKANTU_USE_MPI)
    const StaticCommunicatorMPI & comm_mpi =
      dynamic_cast<const StaticCommunicatorMPI &>(comm.getRealStaticCommunicator());
    if(mumps_data.comm_fortran == MPI_Comm_c2f(comm_mpi.getMPICommunicator()))
#endif
      destroyMumpsData();
  } catch(...) {}


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::initMumpsData(SolverMumpsOptions::ParallelMethod parallel_method) {
  switch(parallel_method) {
  case SolverMumpsOptions::_fully_distributed:
      icntl(18) = 3; //fully distributed
      icntl(28) = 0; //automatic choice

      mumps_data.nz_loc  = matrix->getNbNonZero();
      mumps_data.irn_loc = matrix->getIRN().values;
      mumps_data.jcn_loc = matrix->getJCN().values;
      break;
  case SolverMumpsOptions::_master_slave_distributed:
    if(prank == 0) {
      mumps_data.nz  = matrix->getNbNonZero();
      mumps_data.irn = matrix->getIRN().values;
      mumps_data.jcn = matrix->getJCN().values;
    } else {
      mumps_data.nz  = 0;
      mumps_data.irn = NULL;
      mumps_data.jcn = NULL;

      icntl(18) = 0; //centralized
      icntl(28) = 0; //sequential analysis
    }
    break;
  }
}

/* -------------------------------------------------------------------------- */
void SolverMumps::initialize(SolverOptions & options) {
  AKANTU_DEBUG_IN();

  mumps_data.par = 1;

  if(SolverMumpsOptions * opt = dynamic_cast<SolverMumpsOptions *>(&options)) {
    if(opt->parallel_method == SolverMumpsOptions::_master_slave_distributed) {
      mumps_data.par = 0;
    }
  }

  mumps_data.sym = 2 * (matrix->getSparseMatrixType() == _symmetric);
  prank = communicator.whoAmI();
#ifdef AKANTU_USE_MPI
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<const StaticCommunicatorMPI &>(communicator.getRealStaticCommunicator()).getMPICommunicator());
#endif

  if(AKANTU_DEBUG_TEST(dblTrace)) {
    icntl(1) = 2;
    icntl(2) = 2;
    icntl(3) = 2;
    icntl(4) = 4;
  }

  mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&mumps_data);
  is_mumps_data_initialized = true;

  /* ------------------------------------------------------------------------ */
  UInt size = matrix->getSize();

  if(prank == 0) {
    std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";
    rhs = &(alloc<Real>(sstr_rhs.str(), size, 1, REAL_INIT_VALUE));
  } else {
    rhs = NULL;
  }

  /// No outputs
  icntl(1) = 0;
  icntl(2) = 0;
  icntl(3) = 0;
  icntl(4) = 0;


  mumps_data.nz_alloc = 0;

  if (AKANTU_DEBUG_TEST(dblDump)) icntl(4) = 4;

  mumps_data.n   = size;

  if(AKANTU_DEBUG_TEST(dblDump)) {
    strcpy(mumps_data.write_problem, "mumps_matrix.mtx");
  }

  /* ------------------------------------------------------------------------ */
  // Default Scaling
  icntl(8) = 77;

  icntl(5) = 0; // Assembled matrix

  SolverMumpsOptions * opt = dynamic_cast<SolverMumpsOptions *>(&options);
  if(opt)
    parallel_method = opt->parallel_method;

  initMumpsData(parallel_method);

  mumps_data.job = _smj_analyze; //analyze
  dmumps_c(&mumps_data);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::setRHS(Array<Real> & rhs) {
  if(prank == 0) {
    matrix->getDOFSynchronizer().gather(rhs, 0, this->rhs);
  } else {
    matrix->getDOFSynchronizer().gather(rhs, 0);
  }
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve() {
  AKANTU_DEBUG_IN();

  if(parallel_method == SolverMumpsOptions::_fully_distributed)
    mumps_data.a_loc  = matrix->getA().values;
  else
    if(prank == 0) {
      mumps_data.a  = matrix->getA().values;
    }

  if(prank == 0) {
    mumps_data.rhs = rhs->values;
  }

  /// Default centralized dense second member
  icntl(20) = 0;
  icntl(21) = 0;

  mumps_data.job = _smj_factorize_solve; //solve
  dmumps_c(&mumps_data);

  AKANTU_DEBUG_ASSERT(info(1) != -10, "Singular matrix");
  AKANTU_DEBUG_ASSERT(info(1) == 0,
		      "Error in mumps during solve process, check mumps user guide INFO(1) ="
		      << info(1));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve(Array<Real> & solution) {
  AKANTU_DEBUG_IN();

  solve();

  if(prank == 0) {
    matrix->getDOFSynchronizer().scatter(solution, 0, this->rhs);
  } else {
    matrix->getDOFSynchronizer().scatter(solution, 0);
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
