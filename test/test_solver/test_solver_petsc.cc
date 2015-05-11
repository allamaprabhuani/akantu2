/**
 * @file   test_solver_petsc.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue Dec  2 17:17:34 2014
 *
 * @brief  test the PETSc solver interface
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
/* -------------------------------------------------------------------------- */
#include <cstdlib>
/* -------------------------------------------------------------------------- */
#include "static_communicator.hh"
#include "aka_common.hh"
#include "aka_csr.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "distributed_synchronizer.hh"
#include "petsc_matrix.hh"
#include "solver_petsc.hh"
#include "fe_engine.hh"
#include "dof_synchronizer.hh"

#include "mesh_partition_scotch.hh"

using namespace akantu;
int main(int argc, char *argv[]) {

  initialize(argc, argv);
  const ElementType element_type = _segment_2;
  const GhostType ghost_type = _not_ghost; 
  UInt spatial_dimension = 1;


  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partition it
  Mesh mesh(spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  DistributedSynchronizer * communicator = NULL;
  if(prank == 0) {
    /// creation mesh
    mesh.read("1D_bar.msh");
    MeshPartitionScotch * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }
  

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange,_ek_regular>(mesh, spatial_dimension, "my_fem");

  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  dof_synchronizer.initGlobalDOFEquationNumbers();

  // fill the matrix with 
  UInt nb_element = mesh.getNbElement(element_type);
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element_type);
  UInt nb_dofs_per_element = spatial_dimension * nb_nodes_per_element;
  SparseMatrix K(nb_global_nodes * spatial_dimension, _symmetric);
  K.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  Matrix<Real> element_input(nb_dofs_per_element, nb_dofs_per_element, 0);
  for (UInt i = 0; i < nb_dofs_per_element; ++i) {
    for (UInt j = 0; j < nb_dofs_per_element; ++j) {
      element_input(i, j) = ((i == j) ? 1 : -1); 
    }
  }
  Array<Real> K_e = Array<Real>(nb_element, nb_dofs_per_element * nb_dofs_per_element, "K_e");
  Array<Real>::matrix_iterator K_e_it = K_e.begin(nb_dofs_per_element, nb_dofs_per_element);
  Array<Real>::matrix_iterator K_e_end = K_e.end(nb_dofs_per_element, nb_dofs_per_element);

  for(; K_e_it != K_e_end; ++K_e_it)
    *K_e_it = element_input;

  // assemble the test matrix
  fem->assembleMatrix(K_e, K, spatial_dimension, element_type, ghost_type);

  // apply boundary: block first node
  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();
  Array<bool> boundary = Array<bool>(nb_nodes, spatial_dimension, false);
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (std::abs(position(i, 0)) < Math::getTolerance() )
      boundary(i, 0) = true;
  }

  K.applyBoundary(boundary);

  /// create the PETSc matrix for the solve step
  PETScMatrix petsc_matrix(nb_global_nodes * spatial_dimension, _symmetric);
  petsc_matrix.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  
  /// copy the stiffness matrix into the petsc matrix
  petsc_matrix.add(K, 1);

  // initialize internal forces: they are zero because imposed displacement is zero
  Array<Real> internal_forces(nb_nodes, spatial_dimension, 0.);

  // compute residual: apply nodal force on last node
  Array<Real> residual(nb_nodes, spatial_dimension, 0.);
  
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (std::abs(position(i, 0) - 10) < Math::getTolerance() )
      residual(i, 0) += 2;
  }

  residual -= internal_forces;
  
  /// initialize solver and solution
  Array<Real> solution(nb_nodes, spatial_dimension, 0.);
  SolverPETSc solver(petsc_matrix);
  solver.initialize();
  solver.setOperators();
  solver.setRHS(residual);
  solver.solve(solution);

  /// verify solution
  Math::setTolerance(1e-11);
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (!dof_synchronizer.isPureGhostDOF(i) && !Math::are_float_equal(2 * position(i, 0), solution(i, 0))) {
      std::cout << "The solution is not correct!!!!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  delete communicator;

  finalize();
  
  return EXIT_SUCCESS;

}
