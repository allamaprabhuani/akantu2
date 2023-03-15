/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "fe_engine.hh"
#include "integrator_gauss_igfem.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

void interpolate(const ElementType type);
void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename);

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  /// interpolate test of _igfem_triangle_4
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_4...");
  interpolate(_igfem_triangle_4);

  /// interpolate test of _igfem_triangle_5
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_5...");
  interpolate(_igfem_triangle_5);

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void interpolate(const ElementType type) {

  Int dim = 2;
  std::stringstream mesh_info;
  mesh_info << "mesh_info" << type << ".txt";
  Mesh my_mesh(dim);
  generateIGFEMMesh(type, my_mesh, mesh_info.str());

  Real eps = 3e-13;
  std::cout << "Epsilon : " << eps << std::endl;

  FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem> * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>(
          my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  Array<Real> val_on_quad(nb_quadrature_points, 2, "val_on_quad");

  /// the number of standard nodes in the mesh, i.e. not interface nodes
  UInt nb_standard_nodes = 9;

  /// impose constant value at standard nodes
  for (Int i = 0; i < nb_standard_nodes; ++i) {
    const_val.data()[i * 2 + 0] = 1.;
    const_val.data()[i * 2 + 1] = 2.;
  }

  /// for field to be constant the enriched values need to be zero,
  /// because enrichment is not needed since there is no kink in the
  /// applied field
  for (Int i = nb_standard_nodes; i < const_val.getSize(); ++i) {
    const_val.data()[i * 2 + 0] = 0.;
    const_val.data()[i * 2 + 1] = 0.;
  }

  fem->interpolateOnIntegrationPoints(const_val, val_on_quad, 2, type);

  std::cout << "Interpolation of array : " << const_val << std::endl;
  std::cout << "Gives on quads : " << val_on_quad << std::endl;

  /// interpolate coordinates
  Array<Real> coord_on_quad(nb_quadrature_points, my_mesh.getSpatialDimension(),
                            "coord_on_quad");
  /// create an array with the nodal coordinates that need to be
  /// interpolated. The nodal coordinates of the enriched nodes need
  /// to be set to zero, because they represent the enrichment of the
  /// position field, and the enrichments for this field are all zero!
  /// There is no gap in the mesh!
  Array<Real> igfem_nodes(my_mesh.getNbNodes(), dim);
  const ShapeLagrange<_ek_igfem> & shapes = fem->getShapeFunctions();
  shapes.extractValuesAtStandardNodes(my_mesh.getNodes(), igfem_nodes,
                                      _not_ghost);
  fem->interpolateOnIntegrationPoints(igfem_nodes, coord_on_quad,
                                      my_mesh.getSpatialDimension(), type);

  std::cout << "Interpolations of node coordinates : " << my_mesh.getNodes()
            << std::endl;
  std::cout << "Gives : " << coord_on_quad << std::endl;

  delete fem;
}
