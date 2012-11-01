/**
 * @file   test_cohesive_buildfacets_tetrahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Oct 03 10:20:53 2012
 *
 * @brief  Test for cohesive elements
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
#include <iostream>
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  //  debug::setDebugLevel(dblDump);

  const UInt spatial_dimension = 3;
  const ElementType type = _tetrahedron_10;

  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("tetrahedron.msh", mesh);
  Mesh mesh_facets(spatial_dimension, const_cast<Vector<Real> &>(mesh.getNodes()), "mesh_facets", 1);

  MeshUtils::buildAllFacets(mesh, mesh_facets);

  //  std::cout << mesh << std::endl;

  const ElementType type_facet = mesh.getFacetElementType(type);
  const ElementType type_subfacet = mesh.getFacetElementType(type_facet);
  const ElementType type_subsubfacet = mesh.getFacetElementType(type_subfacet);

  /* ------------------------------------------------------------------------ */
  /* Element to Subelement testing                                            */
  /* ------------------------------------------------------------------------ */

  const Vector<Vector<Element> > & el_to_subel3 = mesh_facets.getElementToSubelement(type_facet);
  const Vector<Vector<Element> > & el_to_subel2 = mesh_facets.getElementToSubelement(type_subfacet);
  const Vector<Vector<Element> > & el_to_subel1 = mesh_facets.getElementToSubelement(type_subsubfacet);

  /// build vectors for comparison
  Vector<Element> tetrahedron(2);
  tetrahedron(0).type = type;
  tetrahedron(0).element = 1;

  tetrahedron(1).type = type;
  tetrahedron(1).element = 11;

  Vector<Element> triangle(8);
  triangle(0).type = type_facet;
  triangle(0).element = 2;

  triangle(1).type = type_facet;
  triangle(1).element = 0;

  triangle(2).type = type_facet;
  triangle(2).element = 4;

  triangle(3).type = type_facet;
  triangle(3).element = 7;

  triangle(4).type = type_facet;
  triangle(4).element = 26;

  triangle(5).type = type_facet;
  triangle(5).element = 24;

  triangle(6).type = type_facet;
  triangle(6).element = 16;

  triangle(7).type = type_facet;
  triangle(7).element = 18;

  Vector<Element> segment(13);
  segment(0).type = type_subfacet;
  segment(0).element = 0;

  segment(1).type = type_subfacet;
  segment(1).element = 1;

  segment(2).type = type_subfacet;
  segment(2).element = 3;

  segment(3).type = type_subfacet;
  segment(3).element = 7;

  segment(4).type = type_subfacet;
  segment(4).element = 9;

  segment(5).type = type_subfacet;
  segment(5).element = 12;

  segment(6).type = type_subfacet;
  segment(6).element = 13;

  segment(7).type = type_subfacet;
  segment(7).element = 16;

  segment(8).type = type_subfacet;
  segment(8).element = 18;

  segment(9).type = type_subfacet;
  segment(9).element = 21;

  segment(10).type = type_subfacet;
  segment(10).element = 27;

  segment(11).type = type_subfacet;
  segment(11).element = 32;

  segment(12).type = type_subfacet;
  segment(12).element = 34;

  /// comparison

  for (UInt i = 0; i < tetrahedron.getSize(); ++i) {
    if (tetrahedron(i).type != el_to_subel3(4)(i).type ||
	tetrahedron(i).element != el_to_subel3(4)(i).element) {
      std::cout << tetrahedron(i).element << " " << el_to_subel3(4)(i).element << std::endl;
      std::cout << "The two tetrahedrons connected to triangle 4 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < triangle.getSize(); ++i) {
    if (triangle(i).type != el_to_subel2(0)(i).type ||
	triangle(i).element != el_to_subel2(0)(i).element) {
      std::cout << "The triangles connected to segment 0 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < segment.getSize(); ++i) {
    if (segment(i).type != el_to_subel1(1)(i).type ||
	segment(i).element != el_to_subel1(1)(i).element) {
      std::cout << "The segments connected to point 1 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }


  /* ------------------------------------------------------------------------ */
  /* Subelement to Element testing                                            */
  /* ------------------------------------------------------------------------ */

  const Vector<Element> & subel_to_el3 = mesh_facets.getSubelementToElement(type);
  const Vector<Element> & subel_to_el2 = mesh_facets.getSubelementToElement(type_facet);
  const Vector<Element> & subel_to_el1 = mesh_facets.getSubelementToElement(type_subfacet);

  /// build vectors for comparison
  Vector<Element> triangle2(mesh.getNbFacetsPerElement(type));
  triangle2(0).type = type_facet;
  triangle2(0).element = 4;

  triangle2(1).type = type_facet;
  triangle2(1).element = 5;

  triangle2(2).type = type_facet;
  triangle2(2).element = 6;

  triangle2(3).type = type_facet;
  triangle2(3).element = 7;

  Vector<Element> segment2(3);
  segment2(0).type = type_subfacet;
  segment2(0).element = 1;

  segment2(1).type = type_subfacet;
  segment2(1).element = 3;

  segment2(2).type = type_subfacet;
  segment2(2).element = 4;

  Vector<Element> point(2);
  point(0).type = mesh.getFacetElementType(type_subfacet);
  point(0).element = 1;

  point(1).type = mesh.getFacetElementType(type_subfacet);
  point(1).element = 2;

  /// comparison

  for (UInt i = 0; i < triangle2.getSize(); ++i) {
    if (triangle2(i).type != subel_to_el3(1, i).type ||
	triangle2(i).element != subel_to_el3(1, i).element) {
      std::cout << "The triangles connected to tetrahedron 1 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < segment2.getSize(); ++i) {
    if (segment2(i).type != subel_to_el2(1, i).type ||
	segment2(i).element != subel_to_el2(1, i).element) {
      std::cout << "The segments connected to triangle 1 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  for (UInt i = 0; i < point.getSize(); ++i) {
    if (point(i).type != subel_to_el1(1, i).type ||
	point(i).element != subel_to_el1(1, i).element) {
      std::cout << "The points connected to segment 1 are wrong"
		<< std::endl;
      return EXIT_FAILURE;
    }
  }


  finalize();

  std::cout << "OK: test_cohesive_buildfacets was passed!" << std::endl;

  return EXIT_SUCCESS;
}
