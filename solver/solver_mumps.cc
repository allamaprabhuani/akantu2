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
#include <mpi.h>

#include "static_communicator_mpi.hh"
#endif

#include "solver_mumps.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SolverMumps::SolverMumps(SparseMatrix & matrix,
			 const SolverID & id,
			 const MemoryID & memory_id) :
  Solver(matrix, id, memory_id) {
  AKANTU_DEBUG_IN();

  UInt size = matrix.getSize();
  //  UInt nb_degre_of_freedom = matrix.getNbDegreOfFreedom();

  //  std::stringstream sstr; sstr << id << ":sparse_matrix";
  //  matrix = new SparseMatrix(mesh, sparse_matrix_type, nb_degre_of_freedom, sstr_mat.str(), memory_id);

  mumps_data.sym = 2 * (matrix.getSparseMatrixType() == _symmetric);

  communicator = StaticCommunicator::getStaticCommunicator();

#ifdef AKANTU_USE_MPI
  mumps_data.par = 1;
  mumps_data.comm_fortran = MPI_Comm_c2f(dynamic_cast<const StaticCommunicatorMPI *>(communicator)->getMPICommunicator());
#else
  mumps_data.par = 0;
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

  mumps_data.job = _smj_initialize; //initialize
  dmumps_c(&mumps_data);

  mumps_data.nz_alloc = 0;

  /// No outputs
  icntl(1) = 0;
  icntl(2) = 0;
  icntl(3) = 0;
  icntl(4) = 0;

  if(AKANTU_DEBUG_TEST(dblTrace)) {
    icntl(1) = 2;
    icntl(2) = 2;
    icntl(3) = 2;
    icntl(4) = 2;
  }

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
void SolverMumps::initNodesLocation(const Mesh & mesh, UInt nb_degre_of_freedom) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_MPI
  if(communicator->getNbProc() > 1) {
    nb_local_nodes = mesh.getNbNodes();
    Vector<UInt> local_nodes(0,2);

    nodes_type = mesh.getNodesType().values;
    UInt * global_node_id = mesh.getGlobalNodesIds().values;

    UInt local_node_val[2];
    for (UInt n = 0; n < nb_local_nodes; ++n) {
      if(nodes_type[n] != -3) {
	local_node_val[0] = global_node_id[n];
	if(nodes_type[n] == -1 || nodes_type[n] == -2) {
	  local_node_val[1] = 0;
	} else if (nodes_type[n] >= -1) {
	  local_node_val[1] = 1;
	}
	local_nodes.push_back(local_node_val);
      }
    }

    nb_local_nodes = local_nodes.getSize();

    UInt nb_proc = communicator->getNbProc();
    if(communicator->whoAmI() == 0) {
      nb_nodes_per_proc[0] = nb_local_nodes;

      communicator->gather(nb_nodes_per_proc, 1);

      for (UInt p = 0; p < nb_proc; ++p) {
	UInt * buffer;
	if(p == 0) buffer = local_nodes.values;
	else {
	  buffer = new UInt[nb_nodes_per_proc[p] * 2];
	  communicator->receive(buffer, 2 * nb_nodes_per_proc[p], p, 0);
	}

	solution_position[p]->resize(0);
	nb_nodes_per_proc_rhs[p] = 0;
	for (UInt n = 0; n < nb_nodes_per_proc[p]; ++n) {
	  UInt node = buffer[2 * n] * nb_degre_of_freedom;
	  if(buffer[2*n + 1] == 0) {
	    nb_nodes_per_proc_rhs[p]++;
	    for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	      rhs_position[p]->push_back(node + d);
	    }
	  }
	  for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	    solution_position[p]->push_back(node + d);
	  }
	}

	if(p != 0) delete [] buffer;
      }
    } else {
      communicator->gather(&nb_local_nodes, 1);
      communicator->send(local_nodes.values, 2 * nb_local_nodes, 0, 0);
    }
  }
#endif // AKANTU_USE_MPI

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

#ifdef AKANTU_USE_MPI
  if(communicator->getNbProc() > 1) {
    Vector<Real> local_rhs(0,1) ;

    Real * rhs_val = rhs.values;
    UInt nb_degre_of_freedom = rhs.getNbComponent();
    UInt nb_nodes = rhs.getSize();
    UInt nb_local_nodes = 0;
    for (UInt n = 0; n < nb_nodes; ++n) {
      if(nodes_type[n] == -1 || nodes_type[n] == -2) {
	nb_local_nodes++;
	UInt node = n * nb_degre_of_freedom;
	for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	  local_rhs.push_back(rhs_val[node + d]);
	}
      }
    }

    Int nb_proc = communicator->getNbProc();
    if (communicator->whoAmI() == 0) {
      this->rhs->clear();
      for (Int p = 0; p < nb_proc; ++p) {
	Real * buffer;
	if(p == 0) buffer = local_rhs.values;
	else {
	  buffer = new Real[nb_degre_of_freedom * nb_nodes_per_proc_rhs[p]];
	  communicator->receive(buffer, nb_degre_of_freedom * nb_nodes_per_proc_rhs[p], p, 0);
	}

	Real * buffer_tmp = buffer;

	for (UInt n = 0; n < nb_nodes_per_proc_rhs[p]; ++n) {
	  UInt node = n * nb_degre_of_freedom;
	  for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	    (*this->rhs)((*rhs_position[p])(node + d)) = *(buffer_tmp++);
	  }
	}
	if(p != 0) delete [] buffer;
      }
    } else {
      communicator->send(local_rhs.values, nb_degre_of_freedom * nb_local_nodes, 0, 0);
    }
  } else {
#endif
  AKANTU_DEBUG_ASSERT(rhs.getSize()*rhs.getNbComponent() == this->rhs->getSize(),
		      "Size of rhs (" << rhs.getSize()*rhs.getNbComponent()
		      << ") and this->rhs (" << this->rhs->getSize()
		      << ") do not match.");

  memcpy(this->rhs->values, rhs.values, this->rhs->getSize() * sizeof(Real));
#ifdef AKANTU_USE_MPI
  }
#endif

  // if(communicator->whoAmI() == 0) {
  //   debug::setDebugLevel(dblDump);
  //   std::cout << *(this->rhs) << std::endl;
  //   debug::setDebugLevel(dblInfo);
  // }
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
		      "Error in mumps suring solve process, check mumps user guide INFO(1) ="
		      << info(0));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolverMumps::solve(Vector<Real> & solution) {
  AKANTU_DEBUG_IN();

  solve();

  // if(communicator->whoAmI() == 0) {
  //   debug::setDebugLevel(dblDump);
  //   std::cout << *(this->rhs) << std::endl;
  //   debug::setDebugLevel(dblInfo);
  // }


#ifdef AKANTU_USE_MPI
  if(communicator->getNbProc() > 1) {
    solution.clear();
    UInt nb_degre_of_freedom = solution.getNbComponent();

    Vector<Real> local_rhs(nb_local_nodes * nb_degre_of_freedom,1);

    if (communicator->whoAmI() == 0) {
      Int nb_proc = communicator->getNbProc();
      for (Int p = 0; p < nb_proc; ++p) {
	Real * buffer;
	if(p == 0) buffer = local_rhs.values;
	else buffer = new Real[nb_degre_of_freedom * nb_nodes_per_proc[p]];

	Real * buffer_tmp = buffer;
	for (UInt n = 0; n < nb_nodes_per_proc[p]; ++n) {
	  UInt node = n * nb_degre_of_freedom;
	  for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	    *(buffer_tmp++) = (*rhs)((*solution_position[p])(node + d));
	  }
	}

	if(p != 0) {
	  communicator->send(buffer, nb_degre_of_freedom * nb_nodes_per_proc[p], p, 0);
	  delete [] buffer;
	}
      }
    } else {
      communicator->receive(local_rhs.values, nb_degre_of_freedom * nb_local_nodes, 0, 0);
    }

    UInt nb_local_nodes = solution.getSize();
    Real * local_rhs_val = local_rhs.values;
    for (UInt n = 0; n < nb_local_nodes; ++n) {
      if(nodes_type[n] >= -2) {
	UInt node = n * nb_degre_of_freedom;
	for (UInt d = 0; d < nb_degre_of_freedom; ++d) {
	  solution.values[node + d] = *(local_rhs_val++);
	}
      }
    }
  } else {
#endif
  memcpy(solution.values, rhs->values, rhs->getSize() * sizeof(Real));
#ifdef AKANTU_USE_MPI
 }
#endif

  // debug::setDebugLevel(dblDump);
  // std::cout << solution << std::endl;
  // debug::setDebugLevel(dblInfo);

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
