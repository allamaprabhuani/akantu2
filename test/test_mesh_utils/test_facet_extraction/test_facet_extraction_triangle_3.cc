/**
 * @file   test_facet_extraction.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
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
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(&argc,&argv);
  const ElementType type = _triangle_3;
  int dim = ElementClass<type>::getSpatialDimension();

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("square.msh", mesh);

  MeshUtils::buildFacets(mesh,1,1);

#ifdef AKANTU_USE_IOHELPER
  const ElementType surf_type = ElementClass<type>::getFacetElementType();

  unsigned int nb_nodes = mesh.getNbNodes();
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction");
  dumper.SetConnectivity((int*)mesh.getConnectivity(type).values,
   			 TRIANGLE1, mesh.getNbElement(type), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  DumperParaview dumper_facet;
  dumper_facet.SetMode(TEXT);

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_boundary");
  dumper_facet.SetConnectivity((int*)mesh.getConnectivity(surf_type).values,
			       LINE1, mesh.getNbElement(surf_type), C_MODE);
  dumper_facet.SetPrefix("paraview/");
  dumper_facet.Init();
  dumper_facet.Dump();

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_internal");
  dumper_facet.SetConnectivity((int*)mesh.getInternalFacetsMesh().getConnectivity(surf_type).values,
  			       LINE1, mesh.getInternalFacetsMesh().getNbElement(surf_type), C_MODE);
  dumper_facet.Init();
  dumper_facet.Dump();

#endif //AKANTU_USE_IOHELPER
  akantu::finalize();
  return EXIT_SUCCESS;
}
