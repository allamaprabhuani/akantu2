/**
 * @file   test_interpolate_tetrahedron_4.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jul 19 10:55:49 2010
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fem.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"


/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  UInt dim = 3;
  ElementType type = _tetrahedron_10;
  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);
  mesh_io.read("cube2.msh", my_mesh);
  FEM *fem = new FEM(my_mesh, dim, "my_fem");

  debug::setDebugLevel(dblDump);
  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  StaticMemory * st_mem = StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  Vector<Real> const_val(fem->getMesh().getNbNodes(), dim, "const_val");
  Vector<Real> val_on_quad(0, dim, "val_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    for (UInt j = 0; j < dim; ++j) {
      const_val.values[i * dim + j] = j;
    }
  }

  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, dim, type);
  std::ofstream my_file("out.txt");
  my_file << const_val << std::endl;
  my_file << val_on_quad << std::endl;

  // interpolate coordinates
  Vector<Real> coord_on_quad(0, my_mesh.getSpatialDimension(), "coord_on_quad");
  fem->interpolateOnQuadraturePoints(my_mesh.getNodes(), coord_on_quad, my_mesh.getSpatialDimension(), type);
  my_file << my_mesh.getNodes() << std::endl;
  my_file << coord_on_quad << std::endl;

  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
