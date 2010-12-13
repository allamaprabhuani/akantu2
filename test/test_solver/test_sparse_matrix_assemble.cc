/**
 * @file   test_sparse_matrix_assemble.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Nov  5 11:13:33 2010
 *
 * @brief  test the assembling method of the SparseMatrix class
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


  akantu::SparseMatrix sparse_matrix(mesh, akantu::_symmetric, spatial_dimension);

  sparse_matrix.buildProfile();

  const akantu::Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  akantu::Mesh::ConnectivityTypeList::const_iterator it;

  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(mesh.getSpatialDimension(*it) != spatial_dimension) continue;
    akantu::UInt nb_element           = mesh.getNbElement(*it);
    akantu::UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*it);
    akantu::Element element(*it);

    akantu::UInt m = nb_nodes_per_element * spatial_dimension;
    akantu::Vector<akantu::Real> local_mat(m, m, 1, "local_mat");

    for(akantu::UInt e = 0; e < nb_element; ++e) {
      element.element = e;
      sparse_matrix.addToMatrix(local_mat, element);
    }
  }

  sparse_matrix.saveMatrix("matrix.mtx");

  akantu::finalize();

  return EXIT_SUCCESS;
}
