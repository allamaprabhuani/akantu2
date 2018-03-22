/**
 * @file   test_mesh_communicator.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Wed Mar 21 2018
 *
 * @brief A Documented file.
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
#include "mesh.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#if defined(AKANTU_USE_MPI)
#  include <mpi.h>
#endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

TEST(TestMeshCommunicator, Default) {
  Mesh mesh(3);

  const Communicator & comm = mesh.getCommunicator();
  ASSERT_EQ(comm.whoAmI(), 0);
  ASSERT_EQ(comm.getNbProc(), 1);
}

TEST(TestMeshCommunicator, Distributed) {
  Mesh mesh(3);
  mesh.distribute();

  int prank = 0, psize = 1;
#if defined(AKANTU_USE_MPI)
  MPI_Comm_size(MPI_COMM_WORLD, &psize);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
#endif
  const Communicator & comm = mesh.getCommunicator();
  ASSERT_EQ(comm.whoAmI(), prank);
  ASSERT_EQ(comm.getNbProc(), psize);
}
