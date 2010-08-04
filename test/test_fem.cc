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
#include "aka_common.hh"
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

  akantu::Vector<akantu::Real> const_val(fem->getNbNodes(), 2, "const_val");
  akantu::Vector<akantu::Real> val_on_quad(0, 2, "val_on_quad");
  akantu::Vector<akantu::Real> grad_on_quad(0, 2, "grad_on_quad");
  akantu::Vector<akantu::Real> int_val_on_elem(0, 2, "int_val_on_elem");
  akantu::Vector<akantu::Real> val_on_nodes_per_elem(fem->getNbElement(akantu::_line_1), 2 * 2,"val_on_nodes_per_elem");
  akantu::Vector<akantu::Real> int_val_on_nodes_per_elem(0, 2 * 2,"int_val_on_nodes_per_elem");
  akantu::Vector<akantu::Real> assemble_val_on_nodes(0, 2,"assemble_val_on_nodes");

  const akantu::Vector<akantu::Real> & shapes = fem->getShapes(akantu::_line_1);

  for (akantu::UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, akantu::_line_1);
  std::cout << const_val << std::endl;
  std::cout << val_on_quad << std::endl;

  fem->gradientOnQuadraturePoints(const_val, grad_on_quad, 2, akantu::_line_1);
  std::cout << grad_on_quad << std::endl;

  fem->integrate(val_on_quad, int_val_on_elem, 2, akantu::_line_1);
  std::cout << int_val_on_elem << std::endl;

  for (akantu::UInt el = 0; el < shapes.getSize(); ++el) {
    val_on_nodes_per_elem.values[el * 4 + 0] = val_on_quad.values[el * 2 + 0] * shapes.values[el * 2 + 0];
    val_on_nodes_per_elem.values[el * 4 + 1] = val_on_quad.values[el * 2 + 1] * shapes.values[el * 2 + 0];
    val_on_nodes_per_elem.values[el * 4 + 2] = val_on_quad.values[el * 2 + 0] * shapes.values[el * 2 + 1];
    val_on_nodes_per_elem.values[el * 4 + 3] = val_on_quad.values[el * 2 + 1] * shapes.values[el * 2 + 1];

  }
  std::cout << val_on_nodes_per_elem << std::endl;

  fem->integrate(val_on_nodes_per_elem, int_val_on_nodes_per_elem, 4, akantu::_line_1);
  std::cout << int_val_on_nodes_per_elem << std::endl;

  fem->assembleVector(int_val_on_nodes_per_elem, assemble_val_on_nodes, 2, akantu::_line_1);
  std::cout << assemble_val_on_nodes << std::endl;

  delete fem;

  akantu::finalize();

  return EXIT_SUCCESS;
}
