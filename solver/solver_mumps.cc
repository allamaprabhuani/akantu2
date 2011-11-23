/**
 * @file   solver_mumps.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 17:32:27 2010
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
  Solver(matrix, id, memory_id), rhs_is_local(true) {
  AKANTU_DEBUG_IN();

  UInt size = matrix.getSize();
  //  UInt nb_degre_of_freedom = matrix.getNbDegreOfFreedom();

  //  std::stringstream sstr; sstr << id << ":sparse_matrix";
  //  matrix = new SparseMatrix(mesh, sparse_matrix_type, nb_degre_of_freedom, sstr_mat.str(), memory_id);

  mumps_data.sym = 2 * (matrix.getSparseMatrixType() == _symmetric);

  communicator = StaticCommunicator::getStaticCommunicator();

  mumps_data.par = 1;
#ifdef AKANTU_USE_MPI
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<const StaticCommunicatorMPI &>(communicator->getRealStaticCommunicator()).getMPICommunicator());
#endif

  if(communicator->whoAmI() == 0) {
    std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";
    rhs = &(alloc<Real>(sstr_rhs.str(), size, 1, REAL_INIT_VALUE));

#ifdef AKANTU_USE_MPI
    UInt nb_proc = communicator->getNbProc();
    nb_nodes_per_proc     = new UInt[nb_proc];
    nb_nodes_per_proc_rhs = new UInt[nb_proc];

    rhs_position      = new Vector<UInt>*[nb_proc];
    solution_position = new Vector<UInt>*[nb_proc];
    for (Int p = 0; p < communicator->getNbProc(); ++p) {
      rhs_position[p]      = new Vector<UInt>(0,1);
      solution_position[p] = new Vector<UInt>(0,1);
    }
#endif
  } else {
    rhs = NULL;
    nb_nodes_per_proc = NULL;
  }

  if(AKANTU_DEBUG_TEST(dblTrace)) {
    icntl(1) = 2;
    icntl(2) = 2;
    icntl(3) = 2;
    icntl(4) = 4;
  }

  mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&mumps_data);


  /// No outputs
  icntl(1) = 0;
  icntl(2) = 0;
  icntl(3) = 0;
  icntl(4) = 0;

  mumps_data.nz_alloc = 0;

  if (AKANTU_DEBUG_TEST(dblDump)) icntl(4) = 4;
  // else if (debug::getDebugLevel() >= dblAccessory)
  //   icntl(4) = 0;
  // else if (debug::getDebugLevel() >= dblSecondary)
  //   icntl(4) = 0;
  // else if (debug::getDebugLevel() >= dblCritical)
  //   icntl(4) = 0;


  mumps_data.n   = size;
  strcpy(mumps_data.write_problem, "mumps_matrix.mtx");


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::~SolverMumps() {
  AKANTU_DEBUG_IN();

  mumps_data.job = _smj_destroy; // destroy
  dmumps_c(&mumps_data);

#ifdef AKANTU_USE_MPI
  if(communicator->whoAmI() == 0) {
    delete [] nb_nodes_per_proc;
    delete [] nb_nodes_per_proc_rhs;

    for (Int p = 0; p < communicator->getNbProc(); ++p) {
      delete rhs_position[p];
      delete solution_position[p];
    }

    delete [] rhs_position;
    delete [] solution_position;
  }
#endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::initialize() {
  AKANTU_DEBUG_IN();

  /// Default Scaling
  icntl(8) = 77;

  icntl(5) = 0; // Assembled matrix

#ifdef AKANTU_USE_MPI
  icntl(18) = 3; //fully distributed
  icntl(28) = 0; //parallel analysis

  mumps_data.nz_loc  = matrix->getNbNonZero();
  mumps_data.irn_loc = matrix->getIRN().values;
  mumps_data.jcn_loc = matrix->getJCN().values;

// #ifdef AKANTU_USE_PTSCOTCH
//   mumps_data.nz_loc  = matrix->getNbNonZero();
//   mumps_data.irn_loc = matrix->getIRN().values;
//   mumps_data.jcn_loc = matrix->getJCN().values;

//   icntl(18) = 3; //fully distributed
//   icntl(28) = 2; //parallel analysis
// #else //  AKANTU_USE_PTSCOTCH
//   icntl(18) = 2; // distributed, sequential analysis
//   icntl(28) = 1; //sequential analysis

  // if (communicator->whoAmI() == 0) {
  //   Int * nb_non_zero_loc = new Int[communicator->getNbProc()];
  //   nb_non_zero_loc[0] = matrix->getNbNonZero();

  //   communicator->gather(nb_non_zero_loc, 1, 0);

  //   UInt nb_non_zero = 0;
  //   for (Int i = 0; i < communicator->getNbProc(); ++i) {
  //     nb_non_zero += nb_non_zero_loc[i];
  //   }

  //   Int * irn = new Int[nb_non_zero];
  //   Int * jcn = new Int[nb_non_zero];

  //   memcpy(irn, matrix->getIRN().values, *nb_non_zero_loc * sizeof(int));
  //   memcpy(jcn, matrix->getJCN().values, *nb_non_zero_loc * sizeof(int));

  //   communicator->gatherv(irn, nb_non_zero_loc, 0);
  //   communicator->gatherv(jcn, nb_non_zero_loc, 0);

  //   mumps_data.nz  = nb_non_zero;
  //   mumps_data.irn = irn;
  //   mumps_data.jcn = jcn;
  // } else {
  //   Int nb_non_zero_loc = matrix->getNbNonZero();
  //   communicator->gather(&nb_non_zero_loc, 1, 0);
  //   communicator->gatherv(matrix->getIRN().values, &nb_non_zero_loc, 0);
  //   communicator->gatherv(matrix->getJCN().values, &nb_non_zero_loc, 0);
  // }
// #endif // AKANTU_USE_PTSCOTCH
#else //AKANTU_USE_MPI
  mumps_data.nz  = matrix->getNbNonZero();
  mumps_data.irn = matrix->getIRN().values;
  mumps_data.jcn = matrix->getJCN().values;

  icntl(18) = 0; //centralized
  icntl(28) = 0; //sequential analysis
#endif //AKANTU_USE_MPI

  mumps_data.job = _smj_analyze; //analyse
  dmumps_c(&mumps_data);

// #ifdef AKANTU_USE_MPI
// #ifdef AKANTU_USE_PTSCOTCH
//   delete [] mumps_data.irn;
//   delete [] mumps_data.jcn;
// #endif
// #endif

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::setRHS(Vector<Real> & rhs) {
  if(communicator->whoAmI() == 0) {
    matrix->getDOFSynchronizer().gather(rhs, 0, this->rhs);
  } else {
    matrix->getDOFSynchronizer().gather(rhs, 0);
  }
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve() {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_MPI
  mumps_data.a_loc  = matrix->getA().values;
#else
  mumps_data.a  = matrix->getA().values;
#endif

  if(communicator->whoAmI() == 0) {
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
void SolverMumps::solve(Vector<Real> & solution) {
  AKANTU_DEBUG_IN();

  solve();

  if(communicator->whoAmI() == 0) {
    matrix->getDOFSynchronizer().scatter(solution, 0, this->rhs);
  } else {
    matrix->getDOFSynchronizer().scatter(solution, 0);
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
