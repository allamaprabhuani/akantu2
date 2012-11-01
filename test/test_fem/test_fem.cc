/**
 * @file   test_fem.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Mon Jul 19 10:55:49 2010
 *
 * @brief  test of the fem class
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
#include <cstdlib>
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
  Mesh my_mesh(1);
  mesh_io.read("line1.msh", my_mesh);
  FEM *fem = new FEM(my_mesh,1,"my_fem");

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  StaticMemory * st_mem = StaticMemory::getStaticMemory();
  std::cout << *st_mem << std::endl;

  Vector<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> val_on_quad(0, 2, "val_on_quad");
  Vector<Real> grad_on_quad(0, 2, "grad_on_quad");
  Vector<Real> int_val_on_elem(0, 2, "int_val_on_elem");
  Vector<Real> val_on_nodes_per_elem(fem->getMesh().getNbElement(_segment_2), 2 * 2,"val_on_nodes_per_elem");
  Vector<Real> int_val_on_nodes_per_elem(0, 2 * 2,"int_val_on_nodes_per_elem");
  Vector<Real> assemble_val_on_nodes(0, 2,"assemble_val_on_nodes");

  const Vector<Real> & shapes = fem->getShapes(_segment_2);

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, _segment_2);
  std::cout << const_val << std::endl;
  std::cout << val_on_quad << std::endl;

  fem->gradientOnQuadraturePoints(const_val, grad_on_quad, 2, _segment_2);
  std::cout << grad_on_quad << std::endl;

  fem->integrate(val_on_quad, int_val_on_elem, 2, _segment_2);
  std::cout << int_val_on_elem << std::endl;

  for (UInt el = 0; el < shapes.getSize(); ++el) {
    val_on_nodes_per_elem.values[el * 4 + 0] = val_on_quad.values[el * 2 + 0] * shapes.values[el * 2 + 0];
    val_on_nodes_per_elem.values[el * 4 + 1] = val_on_quad.values[el * 2 + 1] * shapes.values[el * 2 + 0];
    val_on_nodes_per_elem.values[el * 4 + 2] = val_on_quad.values[el * 2 + 0] * shapes.values[el * 2 + 1];
    val_on_nodes_per_elem.values[el * 4 + 3] = val_on_quad.values[el * 2 + 1] * shapes.values[el * 2 + 1];

  }
  std::cout << val_on_nodes_per_elem << std::endl;

  fem->integrate(val_on_nodes_per_elem, int_val_on_nodes_per_elem, 4, _segment_2);
  std::cout << int_val_on_nodes_per_elem << std::endl;

  fem->assembleVector(int_val_on_nodes_per_elem, assemble_val_on_nodes, 2, _segment_2);
  std::cout << assemble_val_on_nodes << std::endl;

  delete fem;

  finalize();

  return EXIT_SUCCESS;
}
