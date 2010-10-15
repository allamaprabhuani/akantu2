/**
 * @file   static_communicator_mpi.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Sep 14 11:37:10 2010
 *
 * @brief  
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "static_communicator_mpi.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

MPI_Op StaticCommunicatorMPI::synchronizer_operation_to_mpi_op[_so_null + 1] = {
  MPI_SUM,
  MPI_MIN,
  MPI_MAX,
  MPI_OP_NULL
};

__END_AKANTU__
