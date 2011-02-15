/**
 * @file   test_gradient_XXXX.cc
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
  mesh_io.read("FILE.msh", my_mesh);
  FEMTemplate<IntegratorGauss,ShapeLagrange> fem(my_mesh, dim, "my_fem");

  debug::setDebugLevel(dblDump);
  fem.initShapeFunctions();

  std::cout << fem << std::endl;

  Vector<Real> const_val(fem.getMesh().getNbNodes(), 2, "const_val");
  Vector<Real> grad_on_quad(0, 2 * dim, "grad_on_quad");

  for (UInt i = 0; i < const_val.getSize(); ++i) {
    const_val.values[i * 2 + 0] = 1.;
    const_val.values[i * 2 + 1] = 2.;
  }

  fem.gradientOnQuadraturePoints(const_val, grad_on_quad, 2, type);
  std::ofstream my_file("out.txt");
  my_file << const_val << std::endl;
  my_file << grad_on_quad << std::endl;

  // compute gradient of coordinates
  Vector<Real> grad_coord_on_quad(0, dim * dim, "grad_coord_on_quad");
  fem.gradientOnQuadraturePoints(my_mesh.getNodes(), grad_coord_on_quad, my_mesh.getSpatialDimension(), type);
  my_file << my_mesh.getNodes() << std::endl;
  my_file << grad_coord_on_quad << std::endl;

  UInt nb_quads = my_mesh.getNbElement(type) * FEM::getNbQuadraturePoints(type);
  Real eps = 25 * std::numeric_limits<Real>::epsilon();
  std::cout << "Epsilon : " << eps << std::endl;
  for (UInt q = 0; q < nb_quads; ++q) {
    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < dim; ++j) {
	if(!(fabs(grad_coord_on_quad.values[q*dim*dim+ i*dim + j] - (i == j)) <= eps)) {
	  std::cerr << "Error on the quad point " << q << std::endl;
	  for (UInt oi = 0; oi < dim; ++oi) {
	    for (UInt oj = 0; oj < dim; ++oj) {
	      std::cout << fabs(grad_coord_on_quad.values[q*dim*dim + i*dim + j] - (i == j)) << " ";
	    }
	    std::cout << std::endl;
	  }
	  exit(EXIT_FAILURE);
	}
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}
