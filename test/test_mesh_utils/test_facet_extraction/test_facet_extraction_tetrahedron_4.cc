/**
 * @file   test_facet_extraction_tetra1.cc
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
#  include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  int dim = 3;

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);
  Mesh mesh_facets(dim, const_cast<Vector<Real> &>(mesh.getNodes()), "mesh_facets", 1);

  MeshUtils::buildAllFacets(mesh, mesh_facets);

#ifdef AKANTU_USE_IOHELPER
  unsigned int nb_nodes = mesh.getNbNodes();
  iohelper::DumperParaview dumper;
  dumper.SetMode(iohelper::TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_tetrahedron_4).values,
			 iohelper::TETRA1, mesh.getNbElement(_tetrahedron_4), iohelper::C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  iohelper::DumperParaview dumper_facet;
  dumper_facet.SetMode(iohelper::TEXT);

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_boundary");
  dumper_facet.SetConnectivity((int*)mesh.getConnectivity(_triangle_3).values,
			 iohelper::TRIANGLE1, mesh.getNbElement(_triangle_3), iohelper::C_MODE);
  dumper_facet.SetPrefix("paraview/");
  dumper_facet.Init();
  dumper_facet.Dump();

  dumper_facet.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-facet-extraction_internal");
  dumper_facet.SetConnectivity((int*)mesh_facets.getConnectivity(_triangle_3).values,
  			       iohelper::TRIANGLE1, mesh_facets.getNbElement(_triangle_3), iohelper::C_MODE);
  dumper_facet.Init();
  dumper_facet.Dump();


#endif //AKANTU_USE_IOHELPER
  akantu::finalize();
  return EXIT_SUCCESS;
}
