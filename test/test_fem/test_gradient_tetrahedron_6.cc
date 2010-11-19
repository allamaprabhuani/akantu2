/**
 * @file   test_gradient_tetrahedron_6.cc
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
  ElementType type = _tetrahedron_6;
  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);
  mesh_io.read("cube1.msh", my_mesh);
  FEM *fem = new FEM(my_mesh, dim, "my_fem");

  debug::setDebugLevel(dblDump);
  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  StaticMemory * st_mem = StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  Vector<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> grad_on_quad(0, 2 * dim, "grad_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem->gradientOnQuadraturePoints(const_val, grad_on_quad, 2, type);
  std::ofstream my_file("out.txt");
  my_file << const_val << std::endl;
  my_file << grad_on_quad << std::endl;

  // compute gradient of coordinates
  Vector<Real> grad_coord_on_quad(0, 9, "grad_coord_on_quad");
  fem->gradientOnQuadraturePoints(my_mesh.getNodes(), grad_coord_on_quad, my_mesh.getSpatialDimension(), type);
  my_file << my_mesh.getNodes() << std::endl;
  my_file << grad_coord_on_quad << std::endl;

  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
