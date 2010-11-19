/**
 * @file   test_interpolate_segment_3.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Sun Oct  3 16:53:59 2010
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

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
  ElementType type = _segment_3;
  UInt dim = 1;

  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);

  mesh_io.read("line2.msh", my_mesh);

  FEM *fem = new FEM(my_mesh, dim, "my_fem");

  //UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);

  debug::setDebugLevel(dblDump);
  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  StaticMemory * st_mem = StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  Vector<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> val_on_quad(0, 2 , "val_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type);
  std::ofstream my_file("out.txt");
  my_file << const_val << std::endl;
  my_file << val_on_quad << std::endl;

  // interpolate coordinates
  Vector<Real> coord_on_quad(0, my_mesh.getSpatialDimension() , "coord_on_quad");

  fem->interpolateOnQuadraturePoints(my_mesh.getNodes(),
				     coord_on_quad,
				     my_mesh.getSpatialDimension(),
				     type);
  my_file << my_mesh.getNodes() << std::endl;
  my_file << coord_on_quad << std::endl;

  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
