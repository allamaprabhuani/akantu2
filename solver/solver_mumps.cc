/**
 * @file   solver_mumps.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 17:32:27 2010
 *
 * @brief  implem of SolverMumps class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
#include <mpi.h>

#include "static_communicator_mpi.hh"
#endif

#include "solver_mumps.hh"

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
SolverMumps::SolverMumps(const Mesh & mesh,
			 const SparseMatrixType & sparse_matrix_type,
			 UInt nb_degre_of_freedom,
			 const SolverID & id,
			 const MemoryID & memory_id) :
  Solver(mesh, sparse_matrix_type, nb_degre_of_freedom, id, memory_id) {
  AKANTU_DEBUG_IN();

  std::stringstream sstr_mat; sstr_mat << id << ":sparse_matrix";
  matrix = new SparseMatrix(mesh, sparse_matrix_type, nb_degre_of_freedom, sstr_mat.str(), memory_id);

  UInt nb_nodes = mesh.getNbNodes();

  std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";

  rhs = &(alloc<Real>(sstr_rhs.str(), nb_nodes * nb_degre_of_freedom, 1, REAL_INIT_VALUE));

  mumps_data.sym = 2 * (sparse_matrix_type == _symmetric);

#ifdef AKANTU_USE_MPI
  mumps_data.par = 1;
  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<StaticCommunicatorMPI *>(comm)->getMPICommunicator());
#else
  mumps_data.par = 0;
#endif

  icntl(1) = 2;
  icntl(2) = 2;
  icntl(3) = 2;

  icntl(4) = 1;
  if (debug::getDebugLevel() >= 10)
    icntl(4) = 4;
  else if (debug::getDebugLevel() >= 5)
    icntl(4) = 3;
  else if (debug::getDebugLevel() >= 3)
    icntl(4) = 2;

  mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&mumps_data);

  mumps_data.n   = nb_nodes * nb_degre_of_freedom;
  strcpy(mumps_data.write_problem, "mumps_matrix.mtx");

  // mumps_data.nz  = 0;
  // mumps_data.irn = NULL;
  // mumps_data.jcn = NULL;
  // mumps_data.a   = NULL;
  // mumps_data.nz_loc  = 0;
  // mumps_data.irn_loc = NULL;
  // mumps_data.jcn_loc = NULL;
  // mumps_data.a_loc   = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::SolverMumps(SparseMatrix & matrix,
			 const SolverID & id,
			 const MemoryID & memory_id) :
  Solver(matrix, id, memory_id) {
  AKANTU_DEBUG_IN();

  //  std::stringstream sstr; sstr << id << ":sparse_matrix";
  //  matrix = new SparseMatrix(mesh, sparse_matrix_type, nb_degre_of_freedom, sstr_mat.str(), memory_id);

  UInt size = matrix.getSize();
  UInt nb_degre_of_freedom = matrix.getNbDegreOfFreedom();

  std::stringstream sstr_rhs; sstr_rhs << id << ":rhs";


  mumps_data.sym = 2 * (matrix.getSparseMatrixType() == _symmetric);

  communicator = StaticCommunicator::getStaticCommunicator();

#ifdef AKANTU_USE_MPI
  mumps_data.par = 1;
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<const StaticCommunicatorMPI *>(communicator)->getMPICommunicator());
#else
  mumps_data.par = 0;
#endif

  if(communicator->whoAmI() == 0) {
    rhs = &(alloc<Real>(sstr_rhs.str(), size * nb_degre_of_freedom, 1, REAL_INIT_VALUE));
  } else {
    rhs = NULL;
  }

  icntl(1) = 2;
  icntl(2) = 2;
  icntl(3) = 2;

  icntl(4) = 1;
  if (debug::getDebugLevel() >= 10)
    icntl(4) = 4;
  else if (debug::getDebugLevel() >= 5)
    icntl(4) = 3;
  else if (debug::getDebugLevel() >= 3)
    icntl(4) = 2;

  mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&mumps_data);

  mumps_data.n   = size * nb_degre_of_freedom;
  strcpy(mumps_data.write_problem, "mumps_matrix.mtx");

  // mumps_data.nz  = 0;
  // mumps_data.irn = NULL;
  // mumps_data.jcn = NULL;
  // mumps_data.a   = NULL;
  // mumps_data.nz_loc  = 0;
  // mumps_data.irn_loc = NULL;
  // mumps_data.jcn_loc = NULL;
  // mumps_data.a_loc   = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::~SolverMumps() {
  AKANTU_DEBUG_IN();

  delete matrix;

  mumps_data.job = _smj_destroy; // destroy
  dmumps_c(&mumps_data);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolverMumps::initialize() {
  AKANTU_DEBUG_IN();

  //  matrix->buildProfile();
  icntl(5) = 0; // Assembled matrix

#ifdef AKANTU_USE_MPI
  mumps_data.nz_loc  = matrix->getNbNonZero();
  mumps_data.irn_loc = matrix->getIRN().values;
  mumps_data.jcn_loc = matrix->getJCN().values;

  icntl(18) = 3; //fully distributed
  icntl(28) = 0; //parallel analysis

// #ifdef AKANTU_USE_PTSCOTCH
//   icntl(18) = 3; //fully distributed
//   icntl(28) = 2; //parallel analysis
// #else //  AKANTU_USE_PTSCOTCH
//   icntl(18) = 2; // distributed, sequential analysis
//   icntl(28) = 1; //sequential analysis

//   if (communicator->whoAmI() == 0) {
//     Int * nb_non_zero_loc = new Int[communicator->getNbProc()];
//     nb_non_zero_loc[0] = matrix->getNbNonZero();

//     communicator->gather(nb_non_zero_loc, 1, 0);

//     UInt nb_non_zero = 0;
//     for (Int i = 0; i < communicator->getNbProc(); ++i) {
//       nb_non_zero += nb_non_zero_loc[i];
//     }

//     Int * irn = new Int[nb_non_zero];
//     Int * jcn = new Int[nb_non_zero];

//     memcpy(irn, matrix->getIRN().values, *nb_non_zero_loc * sizeof(int));
//     memcpy(jcn, matrix->getJCN().values, *nb_non_zero_loc * sizeof(int));

//     communicator->gatherv(irn, nb_non_zero_loc, 0);
//     communicator->gatherv(jcn, nb_non_zero_loc, 0);

//     mumps_data.nz  = nb_non_zero;
//     mumps_data.irn = irn;
//     mumps_data.jcn = jcn;
//   } else {
//     Int nb_non_zero_loc = matrix->getNbNonZero();
//     communicator->gather(&nb_non_zero_loc, 1, 0);
//     communicator->gatherv(matrix->getIRN().values, &nb_non_zero_loc, 0);
//     communicator->gatherv(matrix->getJCN().values, &nb_non_zero_loc, 0);
//   }
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


#ifdef AKANTU_USE_MPI
#ifdef AKANTU_USE_PTSCOTCH
  delete [] mumps_data.irn;
  delete [] mumps_data.jcn;
#endif
#endif

  AKANTU_DEBUG_OUT();
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
  icntl(20) = 0;
  icntl(21) = 0;

  //  mumps_data.nrhs = 1;
  //  mumps_data.lrhs = mumps_data.n;

  mumps_data.job = _smj_factorize_solve; //solve
  dmumps_c(&mumps_data);

  /// @todo spread the rhs vector form host to slaves

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
