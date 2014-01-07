/**
 * @file   test_facet_extraction_triangle_3.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Fri Sep 03 15:56:15 2010
 *
 * @brief  test of internal facet extraction
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  const ElementType type = _triangle_6;
  int dim = ElementClass<type>::getSpatialDimension();

  Mesh mesh(dim);
  mesh.read("square.msh");
  Mesh mesh_facets(mesh.initMeshFacets("mesh_facets"));

  MeshUtils::buildAllFacets(mesh, mesh_facets);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper1("test-facet-extraction");
  dumper1.registerMesh(mesh, dim);
  dumper1.dump();

  DumperParaview dumper2("test-facet-extraction_boundary");
  dumper2.registerMesh(mesh, dim - 1);
  dumper2.dump();

  DumperParaview dumper3("test-facet-extraction_internal");
  dumper3.registerMesh(mesh_facets, dim);
  dumper3.dump();
#endif //AKANTU_USE_IOHELPER
  akantu::finalize();
  return EXIT_SUCCESS;
}
