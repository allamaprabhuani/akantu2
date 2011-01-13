/**
 * @file   test_gradient_tetrahedron_10.cc
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
  ElementType type = _tetrahedron_10;
  MeshIOMSH mesh_io;
  Mesh my_mesh(dim);
  mesh_io.read("cube2.msh", my_mesh);
  FEM fem(my_mesh, dim, "my_fem");

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
  Vector<Real> grad_coord_on_quad(0, 9, "grad_coord_on_quad");
  fem.gradientOnQuadraturePoints(my_mesh.getNodes(), grad_coord_on_quad, my_mesh.getSpatialDimension(), type);
  my_file << my_mesh.getNodes() << std::endl;
  my_file << grad_coord_on_quad << std::endl;

  UInt nb_quads = my_mesh.getNbElement(type) * FEM::getNbQuadraturePoints(type);
  Real eps = 30 * std::numeric_limits<Real>::epsilon();
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
