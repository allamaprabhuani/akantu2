/**
 * @file   test_integrate_XXXX.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
  ElementType type = TYPE;
  UInt dim = DIM;

  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);

  mesh_io.read("XXXX.msh", my_mesh);

  FEM fem(my_mesh, dim, "my_fem");

  debug::_debug_level = dblDump;
  fem.initShapeFunctions();



  std::cout << fem << std::endl;

  Vector<Real> const_val(fem.getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> val_on_quad(0, 2 , "val_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  //interpolate function on quadrature points
  fem.interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type);
  //integrate function on elements
  akantu::Vector<akantu::Real> int_val_on_elem(0, 2, "int_val_on_elem");
  fem.integrate(val_on_quad, int_val_on_elem, 2, type);

  // get global integration value
  Real value[2] = {0,0};
  std::ofstream my_file("out.txt");
  my_file << val_on_quad << std::endl << int_val_on_elem << std::endl;
  for (UInt i = 0; i < fem.getMesh().getNbElement(type); ++i) {
    value[0] += int_val_on_elem.values[2*i];
    value[1] += int_val_on_elem.values[2*i+1];
  }

  my_file << "integral is " << value[0] << " " << value[1] << std::endl;

  FEM fem_boundary(my_mesh, dim-1, "my_fem_boundary");
  fem_boundary.initShapeFunctions();

  ElementType bound_type = Mesh::getFacetElementType(type);
  UInt nb_boundary_quad  = FEM::getNbQuadraturePoints(bound_type);

  Vector<Real> val_on_bquad(0, nb_boundary_quad, "val_on_quad");
  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
  }

  finalize();

  return EXIT_SUCCESS;
}
