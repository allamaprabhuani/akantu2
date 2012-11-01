/**
 * @file   test_surface_extraction_3d.cc
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Mon Oct 25 11:40:12 2010
 *
 * @brief 3d surface extraction tests
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"

#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  akantu::initialize(argc, argv);

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cubes.msh", mesh);

  MeshUtils::buildFacets(mesh);
  MeshUtils::buildSurfaceID(mesh);

  unsigned int nb_nodes = mesh.getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  dumper.setMode(iohelper::TEXT);

  dumper.setPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction");
  dumper.setConnectivity((int*)mesh.getConnectivity(_tetrahedron_4).values,
   			 iohelper::TETRA1, mesh.getNbElement(_tetrahedron_4), iohelper::C_MODE);
  dumper.setPrefix("paraview/");
  dumper.init();
  dumper.dump();

  iohelper::DumperParaview dumper_surface;
  dumper_surface.setMode(iohelper::TEXT);

  dumper_surface.setPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_boundary");

  dumper_surface.setConnectivity((int *)mesh.getConnectivity(_triangle_3).values,
  			       iohelper::TRIANGLE1, mesh.getNbElement(_triangle_3), iohelper::C_MODE);
  double * surf_id = new double [mesh.getSurfaceID(_triangle_3).getSize()];
  for (UInt i = 0; i < mesh.getSurfaceID(_triangle_3).getSize(); ++i)
    surf_id[i] = (double)mesh.getSurfaceID(_triangle_3).values[i];
  dumper_surface.addElemDataField("surface_id", surf_id, 1, mesh.getNbElement(_triangle_3));
  dumper_surface.setPrefix("paraview/");
  dumper_surface.init();
  dumper_surface.dump();

  delete [] surf_id;
#endif //AKANTU_USE_IOHELPER

  finalize();
  return EXIT_SUCCESS;
}
