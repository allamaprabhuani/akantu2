/**
 * @file   test_integrate_tetrahedron_6.cc
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
  MeshIOMSH mesh_io;
  Mesh my_mesh(3);
  mesh_io.read("cube.msh", my_mesh);
  FEM *fem = new FEM(my_mesh,3,"my_fem");

  debug::_debug_level = dblDump;
  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  StaticMemory * st_mem = StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  Vector<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> val_on_quad(0, 2, "val_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }
  //interpolate function on quadrature points
  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, _tetrahedron_6);
  //integrate function on elements
  akantu::Vector<akantu::Real> int_val_on_elem(0, 2, "int_val_on_elem");
  fem->integrate(val_on_quad,int_val_on_elem,2,_tetrahedron_6);
  // get global integration value
  Real value[2] = {0,0};
  std::ofstream my_file("out.txt");
  my_file << int_val_on_elem << std::endl;
  for (UInt i = 0; i < fem->getMesh().getNbElement(_tetrahedron_6); ++i) {
    value[0] += int_val_on_elem.values[2*i];
    value[1] += int_val_on_elem.values[2*i+1];
  }
  
  my_file << "integral is " << value[0] << " " << value[1] << std::endl;



  //  delete fem;

  //  finalize();

  return EXIT_SUCCESS;
}
