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
  // akantu::UInt n = 10 * comm->getNbProc();

  // akantu::SparseMatrix * sparse_matrix = new akantu::SparseMatrix(n, akantu::_symmetric, 1, "hand");
  // akantu::Solver * solver = new akantu::SolverMumps(*sparse_matrix);

  // akantu::UInt i_start = comm->whoAmI() * 10;
  // for(akantu::UInt i = i_start; i < i_start + 10; ++i) {
  //   sparse_matrix->addToProfile(i,i);
  //   sparse_matrix->addToMatrix(i, i, 1.);
  // }

  // sparse_matrix->addToProfile(9, 8);
  // sparse_matrix->addToMatrix(9, 8, -1.);
  // sparse_matrix->addToMatrix(9, 9, -1.);

// akantu::SparseMatrix * sparse_matrix = new akantu::SparseMatrix(4, akantu::_symmetric, 1, "hand");
// akantu::Solver * solver = new akantu::SolverMumps(*sparse_matrix);

// sparse_matrix->addToProfile(0, 0);
// sparse_matrix->addToProfile(1, 0);
// sparse_matrix->addToProfile(2, 0);
// sparse_matrix->addToProfile(3, 0);
// sparse_matrix->addToProfile(1, 1);
// sparse_matrix->addToProfile(2, 1);
// sparse_matrix->addToProfile(3, 1);
// sparse_matrix->addToProfile(2, 2);
// sparse_matrix->addToProfile(3, 2);
// sparse_matrix->addToProfile(3, 3);

// sparse_matrix->addToMatrix(0, 0, 0.432692  );
// sparse_matrix->addToMatrix(0, 1, 0.240385  );
// sparse_matrix->addToMatrix(0, 2, -0.240385 );
// sparse_matrix->addToMatrix(0, 3, -0.0480769);
// sparse_matrix->addToMatrix(1, 1, 0.432692  );
// sparse_matrix->addToMatrix(1, 2, -0.0480769);
// sparse_matrix->addToMatrix(1, 3, -0.240385 );
// sparse_matrix->addToMatrix(2, 2, 0.432692  );
// sparse_matrix->addToMatrix(2, 3, 0.240385  );
// sparse_matrix->addToMatrix(3, 3, 0.432692  );


  // akantu::Vector<akantu::Real> * rhs;
  // if(comm->whoAmI() == 0) {
  //   rhs = new akantu::Vector<akantu::Real>(4, 1);
  //   rhs->clear();
  //   rhs->values[0] = 10000;
  //   rhs->values[1] = 10000;
  //   // for(akantu::UInt i = 0; i < n; ++i) {
  //   //   rhs->values[i] = 1.;
  //   // }

  //   solver->setRHS(*rhs);
  //   rhs->clear();
  // }


  akantu::SparseMatrix * sparse_matrix = new akantu::SparseMatrix(8, akantu::_symmetric, 1, "hand");
  sparse_matrix->addToProfile(1-1, 1-1);
  sparse_matrix->addToProfile(1-1, 2-1);
  sparse_matrix->addToProfile(1-1, 3-1);
  sparse_matrix->addToProfile(1-1, 4-1);
  sparse_matrix->addToProfile(1-1, 5-1);
  sparse_matrix->addToProfile(1-1, 6-1);
  sparse_matrix->addToProfile(1-1, 7-1);
  sparse_matrix->addToProfile(1-1, 8-1);
  sparse_matrix->addToProfile(2-1, 2-1);
  sparse_matrix->addToProfile(2-1, 3-1);
  sparse_matrix->addToProfile(2-1, 4-1);
  sparse_matrix->addToProfile(2-1, 5-1);
  sparse_matrix->addToProfile(2-1, 6-1);
  sparse_matrix->addToProfile(2-1, 7-1);
  sparse_matrix->addToProfile(2-1, 8-1);
  sparse_matrix->addToProfile(3-1, 3-1);
  sparse_matrix->addToProfile(3-1, 4-1);
  sparse_matrix->addToProfile(3-1, 5-1);
  sparse_matrix->addToProfile(3-1, 6-1);
  sparse_matrix->addToProfile(3-1, 7-1);
  sparse_matrix->addToProfile(3-1, 8-1);
  sparse_matrix->addToProfile(4-1, 4-1);
  sparse_matrix->addToProfile(4-1, 5-1);
  sparse_matrix->addToProfile(4-1, 6-1);
  sparse_matrix->addToProfile(4-1, 7-1);
  sparse_matrix->addToProfile(4-1, 8-1);
  sparse_matrix->addToProfile(5-1, 5-1);
  sparse_matrix->addToProfile(5-1, 6-1);
  sparse_matrix->addToProfile(5-1, 7-1);
  sparse_matrix->addToProfile(5-1, 8-1);
  sparse_matrix->addToProfile(6-1, 6-1);
  sparse_matrix->addToProfile(6-1, 7-1);
  sparse_matrix->addToProfile(6-1, 8-1);
  sparse_matrix->addToProfile(7-1, 7-1);
  sparse_matrix->addToProfile(7-1, 8-1);
  sparse_matrix->addToProfile(8-1, 8-1);

  sparse_matrix->addToMatrix(1-1, 1-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(1-1, 2-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 3-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 4-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 5-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 6-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(1-1, 8-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 2-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(2-1, 3-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 4-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 5-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 6-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(2-1, 8-1, 0		  );
  sparse_matrix->addToMatrix(3-1, 3-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(3-1, 4-1, 0		  );
  sparse_matrix->addToMatrix(3-1, 5-1, 0.240384615384615  );
  sparse_matrix->addToMatrix(3-1, 6-1, -0.240384615384615 );
  sparse_matrix->addToMatrix(3-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(3-1, 8-1, -0.0480769230769231);
  sparse_matrix->addToMatrix(4-1, 4-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(4-1, 5-1, 0		  );
  sparse_matrix->addToMatrix(4-1, 6-1, 0		  );
  sparse_matrix->addToMatrix(4-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(4-1, 8-1, 0		  );
  sparse_matrix->addToMatrix(5-1, 5-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(5-1, 6-1, -0.0480769230769231);
  sparse_matrix->addToMatrix(5-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(5-1, 8-1, -0.240384615384615 );
  sparse_matrix->addToMatrix(6-1, 6-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(6-1, 7-1, 0		  );
  sparse_matrix->addToMatrix(6-1, 8-1, 0.240384615384615  );
  sparse_matrix->addToMatrix(7-1, 7-1, 0.432692307692308  );
  sparse_matrix->addToMatrix(7-1, 8-1, 0		  );
  sparse_matrix->addToMatrix(8-1, 8-1, 0.432692307692308  );

  akantu::Solver * solver = new akantu::SolverMumps(*sparse_matrix);
  akantu::Vector<akantu::Real> * rhs;
  if(comm->whoAmI() == 0) {
    rhs = new akantu::Vector<akantu::Real>(8, 1);
    rhs->clear();
    rhs->values[2] = 10000;
    rhs->values[4] = 10000;
    // for(akantu::UInt i = 0; i < n; ++i) {
    //   rhs->values[i] = 1.;
    // }

    solver->setRHS(*rhs);
    rhs->clear();
  }

  //  std::stringstream sstr; sstr << "solver_matrix.mtx" << comm->whoAmI();
  //  sparse_matrix->saveMatrix(sstr.str());

  solver->initialize();

  solver->solve(*rhs);

  if(comm->whoAmI() == 0) {
    akantu::debug::setDebugLevel(akantu::dblDump);
    std::cout << *rhs << std::endl;
    akantu::debug::setDebugLevel(akantu::dblWarning);
  }

  delete solver;

  akantu::finalize();

  return EXIT_SUCCESS;
}
