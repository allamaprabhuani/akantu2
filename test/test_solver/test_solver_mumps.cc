/**
 * @file   test_solver_mumps.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Dec 12 20:20:32 2010
 *
 * @brief  simple test of the mumps solver interface
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "solver_mumps.hh"
#include "static_communicator.hh"


/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  akantu::initialize(&argc, &argv);

  akantu::debug::setDebugLevel(akantu::dblDump);

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

  //  sparse_matrix->saveMatrix("solver_matrix.mtx");
  solver->initialize();
  solver->solve();

  if(comm->whoAmI() == 0) {
    std::cout << solver->getRHS() << std::endl;
  }

  delete solver;

  akantu::finalize();

  return EXIT_SUCCESS;
}
