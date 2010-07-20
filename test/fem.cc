/**
 * @file   fem.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Jul 19 10:55:49 2010
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include <cstdlib>
/* -------------------------------------------------------------------------- */
#include "common.hh"
#include "fem.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
  akantu::FEM *fem = new akantu::FEM(1, "my_fem");
  akantu::MeshIOMSH mesh_io;

  mesh_io.read("line.msh", fem->getMesh());

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  akantu::StaticMemory * st_mem = akantu::StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  akantu::Vector<akantu::Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  akantu::Vector<akantu::Real> const_val_on_quad(0, 2, "const_val_on_quad");

  for (akantu::UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem->interpolateOnQuadraturePoints(const_val, const_val_on_quad, akantu::_line_1);
  std::cout << const_val << std::endl;
  std::cout << const_val_on_quad << std::endl;

  delete fem;

  akantu::finalize();

  return EXIT_SUCCESS;
}
