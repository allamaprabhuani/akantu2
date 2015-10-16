/**
 * @file   sparse_solver_mumps.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  implem of SparseSolverMumps class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "dof_manager_default.hh"
#include "sparse_matrix_aij.hh"

#if defined(AKANTU_USE_MPI)
#include "static_communicator_mpi.hh"
#include "mpi_type_wrapper.hh"
#endif

#include "solver_mumps.hh"
#include "dof_synchronizer.hh"

/* -------------------------------------------------------------------------- */
// static std::ostream & operator <<(std::ostream & stream, const DMUMPS_STRUC_C
// & _this) {
//   stream << "DMUMPS Data [" << std::endl;
//   stream << " + job          : " << _this.job          << std::endl;
//   stream << " + par          : " << _this.par          << std::endl;
//   stream << " + sym          : " << _this.sym          << std::endl;
//   stream << " + comm_fortran : " << _this.comm_fortran << std::endl;
//   stream << " + nz           : " << _this.nz           << std::endl;
//   stream << " + irn          : " << _this.irn          << std::endl;
//   stream << " + jcn          : " << _this.jcn          << std::endl;
//   stream << " + nz_loc       : " << _this.nz_loc       << std::endl;
//   stream << " + irn_loc      : " << _this.irn_loc      << std::endl;
//   stream << " + jcn_loc      : " << _this.jcn_loc      << std::endl;
//   stream << "]";
//   return stream;
// }

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SparseSolverMumps::SparseSolverMumps(DOFManagerDefault & dof_manager,
                                     const ID & matrix_id, const ID & id,
                                     const MemoryID & memory_id)
    : SparseSolver(dof_manager, matrix_id, id, memory_id),
      dof_manager(dof_manager), matrix(dof_manager.getMatrix(matrix_id)),
      rhs(dof_manager.getResidual()), master_rhs_solution(0, 1) {
  AKANTU_DEBUG_IN();

  StaticCommunicator & communicator = StaticCommunicator::getStaticCommunicator();
  this->prank = communicator.whoAmI();

#ifdef AKANTU_USE_MPI
  this->parallel_method = _fully_distributed;
#else  // AKANTU_USE_MPI
  this->parallel_method = _not_parallel;
#endif // AKANTU_USE_MPI

  this->mumps_data.par = 1; // The host is part of computations

  switch (this->parallel_method) {
  case _not_parallel:
    break;
  case _master_slave_distributed:
    this->mumps_data.par = 0; // The host is not part of the computations
  case _fully_distributed:
#ifdef AKANTU_USE_MPI
    const StaticCommunicatorMPI & mpi_st_comm =
        dynamic_cast<const StaticCommunicatorMPI &>(
            communicator.getRealStaticCommunicator());

    this->mumps_data.comm_fortran =
        MPI_Comm_c2f(mpi_st_comm.getMPITypeWrapper().getMPICommunicator());
#else
    AKANTU_DEBUG_ERROR(
        "You cannot use parallel method to solve without activating MPI");
#endif
    break;
  }

  this->mumps_data.sym = 2 * (matrix.getMatrixType() == _symmetric);
  this->prank = communicator.whoAmI();

  this->mumps_data.job = _smj_initialize; // initialize
  dmumps_c(&this->mumps_data);

  /* ------------------------------------------------------------------------ */
  /* ------------------------------------------------------------------------ */
  // Output setup
  if (AKANTU_DEBUG_TEST(dblTrace)) {
    icntl(1) = 6;
    icntl(2) = 2;
    icntl(3) = 2;
    icntl(4) = 4;
  } else {
    /// No outputs
    icntl(1) = 6; // error output
    icntl(2) = 0; // dignostics output
    icntl(3) = 0; // informations
    icntl(4) = 0; // no outputs
  }

  if (AKANTU_DEBUG_TEST(dblDump)) {
    strcpy(this->mumps_data.write_problem, "mumps_matrix.mtx");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SparseSolverMumps::~SparseSolverMumps() {
  AKANTU_DEBUG_IN();

  this->mumps_data.job = _smj_destroy; // destroy
  dmumps_c(&this->mumps_data);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::initMumpsData() {
  // Default Scaling
  icntl(8) = 77;

  // Assembled matrix
  icntl(5) = 0;

  /// Default centralized dense second member
  icntl(20) = 0;
  icntl(21) = 0;

  // automatic choice for analysis
  icntl(28) = 0;

  UInt size = matrix.getSize();

  if (prank == 0) {
    this->master_rhs_solution.resize(size);
  }

  this->mumps_data.nz_alloc = 0;
  this->mumps_data.n = size;

  switch (this->parallel_method) {
  case _fully_distributed:
    icntl(18) = 3; // fully distributed

    this->mumps_data.nz_loc = matrix.getNbNonZero();
    this->mumps_data.irn_loc = matrix.getIRN().storage();
    this->mumps_data.jcn_loc = matrix.getJCN().storage();

    break;
  case _not_parallel:
  case _master_slave_distributed:
    icntl(18) = 0; // centralized

    if (prank == 0) {
      this->mumps_data.nz = matrix.getNbNonZero();
      this->mumps_data.irn = matrix.getIRN().storage();
      this->mumps_data.jcn = matrix.getJCN().storage();
    } else {
      this->mumps_data.nz = 0;
      this->mumps_data.irn = NULL;
      this->mumps_data.jcn = NULL;
    }
    break;
  default:
    AKANTU_DEBUG_ERROR("This case should not happen!!");
  }
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::initialize() {
  AKANTU_DEBUG_IN();

  this->analysis();

  //  icntl(14) = 80;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::analysis() {
  AKANTU_DEBUG_IN();

  initMumpsData();

  this->mumps_data.job = _smj_analyze; // analyze
  dmumps_c(&this->mumps_data);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::factorize() {
  AKANTU_DEBUG_IN();

  this->mumps_data.rhs = this->rhs.storage();
  if (parallel_method == _fully_distributed)
    this->mumps_data.a_loc = this->matrix.getA().storage();
  else {
    if (prank == 0)
      this->mumps_data.a = this->matrix.getA().storage();
  }

  this->mumps_data.job = _smj_factorize; // factorize
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::solve() {
  AKANTU_DEBUG_IN();

  // if (prank == 0) {
  //   // deactivate debug messages
  //   DebugLevel dbl = debug::getDebugLevel();
  //   debug::setDebugLevel(dblError);

  //   matrix.getDOFSynchronizer().gather(this->rhs, 0, this->master_rhs_solution);

  //   // reactivate debug messages
  //   debug::setDebugLevel(dbl);
  // } else {
  //   this->matrix.getDOFSynchronizer().gather(this->rhs, 0);
  // }

  if (prank == 0)
    this->mumps_data.rhs = this->master_rhs_solution.storage();

  this->mumps_data.job = _smj_solve; // solve
  dmumps_c(&this->mumps_data);

  this->printError();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SparseSolverMumps::printError() {
  Int _info_v[2];
  _info_v[0] = info(1);  // to get errors
  _info_v[1] = -info(1); // to get warnings
  StaticCommunicator::getStaticCommunicator().allReduce(_info_v, 2, _so_min);
  _info_v[1] = -_info_v[1];

  if (_info_v[0] < 0) { // < 0 is an error
    switch (_info_v[0]) {
    case -10:
      AKANTU_DEBUG_ERROR("The matrix is singular");
      break;
    case -9: {
      icntl(14) += 10;
      if (icntl(14) != 90) {
        // std::cout << "Dynamic memory increase of 10%" << std::endl;
        AKANTU_DEBUG_WARNING(
            "MUMPS dynamic memory is insufficient it will be increased of 10%");
        this->analysis();
        this->factorize();
        this->solve();
      } else {
        AKANTU_DEBUG_ERROR("The MUMPS workarray is too small INFO(2)="
                           << info(2) << "No further increase possible");
        break;
      }
    }
    default:
      AKANTU_DEBUG_ERROR("Error in mumps during solve process, check mumps "
                         "user guide INFO(1) = "
                         << _info_v[1]);
    }
  } else if (_info_v[1] > 0) {
    AKANTU_DEBUG_WARNING("Warning in mumps during solve process, check mumps "
                         "user guide INFO(1) = "
                         << _info_v[1]);
  }
}

__END_AKANTU__
