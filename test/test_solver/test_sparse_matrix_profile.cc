/**
 * @file   test_sparse_matrix_profile.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Nov  5 11:13:33 2010
 *
 * @brief  test the profile generation of the SparseMatrix class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdlib>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"

#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  akantu::initialize(&argc, &argv);

  akantu::UInt spatial_dimension = 2;
  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);


  akantu::SparseMatrix sparse_matrix(mesh, akantu::_symmetric, 2);

  sparse_matrix.buildProfile();

  sparse_matrix.saveProfile("profile.mtx");

  akantu::finalize();

  return EXIT_SUCCESS;
}
