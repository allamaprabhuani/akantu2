/**
 * @file   static_communicator_mpi.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Sep 14 11:37:10 2010
 *
 * @brief  
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
