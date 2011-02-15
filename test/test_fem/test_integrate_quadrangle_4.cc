/**
 * @file   test_integrate_quadrangle_4.cc
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
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "fem.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"


/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  ElementType type = _quadrangle_4;
  UInt dim = 2;

  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);

  mesh_io.read("square_structured1.msh", my_mesh);

  FEMTemplate<IntegratorGauss,ShapeLagrange> fem(my_mesh, dim, "my_fem");

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

  FEMTemplate<IntegratorGauss,ShapeLagrange> fem_boundary(my_mesh, dim-1, "my_fem_boundary");
  fem_boundary.initShapeFunctions();

  ElementType bound_type = Mesh::getFacetElementType(type);
  UInt nb_boundary_quad  = fem_boundary.getNbQuadraturePoints(bound_type);

  Vector<Real> val_on_bquad(0, nb_boundary_quad, "val_on_quad");
  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
  }

  finalize();

  return EXIT_SUCCESS;
}
