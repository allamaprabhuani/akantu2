/**
 * @file   test_interpolate_XXXX.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fem.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  ElementType type = XXXX;
  UInt dim = DIM;

  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);

  mesh_io.read("XXXX.msh", my_mesh);

  FEM *fem = new FEMTemplate<IntegratorGauss,ShapeLagrange>(my_mesh, dim, "my_fem");

  //UInt nb_quadrature_points = FEM::getNbQuadraturePoints(type);

  debug::setDebugLevel(dblDump);
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

  fem->interpolateOnQuadraturePoints(const_val, val_on_quad, 2, type);
  std::ofstream my_file("out.txt");
  my_file << const_val << std::endl;
  my_file << val_on_quad << std::endl;

  // interpolate coordinates
  Vector<Real> coord_on_quad(0, my_mesh.getSpatialDimension(), "coord_on_quad");

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
