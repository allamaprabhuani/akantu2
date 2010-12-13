/**
 * @file   solver_mumps.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 17 17:32:27 2010
 *
 * @brief  implem of SolverMumps class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 * @section DESCRIPTION
 *
 * @subsection Ctrl_param Control parameters
 *
 * ICNTL(1),
 * ICNTL(2),
 * ICNTL(3) : output streams for error, diagnostics, ans global messages
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
#endif

#include "solver_mumps.hh"
#include "static_communicator_mpi.hh"


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
  mumps_data.n   = nb_nodes * nb_degre_of_freedom;

  mumps_data.nz  = 0;
  mumps_data.irn = NULL;
  mumps_data.jcn = NULL;
  mumps_data.a   = NULL;
  mumps_data.nz_loc  = 0;
  mumps_data.irn_loc = NULL;
  mumps_data.jcn_loc = NULL;
  mumps_data.a_loc   = NULL;


#ifdef AKANTU_USE_MPI
  mumps_data.par = 1;
  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator();
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<StaticCommunicatorMPI *>(comm)->getMPICommunicator());
#else
  mumps_data.par = 0;
#endif

  mumps_data.job = -1; //initialize
  dmumps_c(&mumps_data);

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

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SolverMumps::~SolverMumps() {
  AKANTU_DEBUG_IN();

  delete matrix;

  mumps_data.job = -2; // destroy
  dmumps_c(&mumps_data);

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void SolverMumps::initialize() {
  AKANTU_DEBUG_IN();

  matrix->buildProfile();

  icntl(5) = 0;

#ifdef AKANTU_USE_MPI
  mumps_data.nz_loc  = matrix->getNbNonZero();
  mumps_data.irn_loc = matrix->getIRN().values;
  mumps_data.jcn_loc = matrix->getJCN().values;

  icntl(18) = 3; //distributed
  icntl(28) = 2; //parallel analysis
#else
  mumps_data.nz  = matrix->getNbNonZero();
  mumps_data.irn = matrix->getIRN().values;
  mumps_data.jcn = matrix->getJCN().values;

  icntl(18) = 0; //centralized
  icntl(28) = 0; //sequential analysis
#endif

  mumps_data.job = 2; //analyse
  dmumps_c(&mumps_data);

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

  mumps_data.job = 3; //factorize
  dmumps_c(&mumps_data);

  mumps_data.job = 2; //solve
  dmumps_c(&mumps_data);

  /// @todo spread the rhs vector form host to slaves

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
