/**
 * @file   test_solver_mumps.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Dec 12 20:20:32 2010
 *
 * @brief  simple test of the mumps solver interface
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
#include "solver_mumps.hh"
#include "static_communicator.hh"


/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::initialize(&argc, &argv);

  //  akantu::debug::setDebugLevel(akantu::dblDump);

  akantu::StaticCommunicator * comm = akantu::StaticCommunicator::getStaticCommunicator();
  akantu::UInt n = 10 * comm->getNbProc();

  akantu::SparseMatrix * sparse_matrix = new akantu::SparseMatrix(n, akantu::_symmetric, 1, "hand");
  akantu::Solver * solver = new akantu::SolverMumps(*sparse_matrix);

  akantu::UInt i_start = comm->whoAmI() * 10;
  for(akantu::UInt i = i_start; i < i_start + 10; ++i) {
    sparse_matrix->addToProfile(i,i);
    sparse_matrix->addToMatrix(i, i, 1./(i+1));
  }

  if(comm->whoAmI() == 0)
    for(akantu::UInt i = 0; i < n; ++i) {
      solver->getRHS().values[i] = 1.;
    }


  std::stringstream sstr; sstr << "solver_matrix.mtx" << comm->whoAmI();
  sparse_matrix->saveMatrix(sstr.str());

  std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB " << comm->whoAmI() << std::endl;


  solver->initialize();
  solver->solve();

  if(comm->whoAmI() == 0) {
    akantu::debug::setDebugLevel(akantu::dblDump);
    std::cout << solver->getRHS() << std::endl;
    akantu::debug::setDebugLevel(akantu::dblWarning);
  }

  delete solver;

  akantu::finalize();

  return EXIT_SUCCESS;
}
