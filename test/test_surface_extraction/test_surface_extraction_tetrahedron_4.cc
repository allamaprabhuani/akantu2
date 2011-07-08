/**
 * @file   test_surface_extraction_3d.cc
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Mon Oct 25 11:40:12 2010
 *
 * @brief  
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
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  akantu::initialize(&argc, &argv);

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cubes.msh", mesh);
  
  MeshUtils::buildFacets(mesh,1,0);
  MeshUtils::buildSurfaceID(mesh);

  unsigned int nb_nodes = mesh.getNbNodes();
#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_tetrahedron_4).values,
   			 TETRA1, mesh.getNbElement(_tetrahedron_4), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

  DumperParaview dumper_surface;
  dumper_surface.SetMode(TEXT);

  dumper_surface.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-surface-extraction_boundary");
  
  dumper_surface.SetConnectivity((int *)mesh.getConnectivity(_triangle_3).values,
  			       TRIANGLE1, mesh.getNbElement(_triangle_3), C_MODE);
  double * surf_id = new double [mesh.getSurfaceID(_triangle_3).getSize()];
  for (UInt i = 0; i < mesh.getSurfaceID(_triangle_3).getSize(); ++i)
    surf_id[i] = (double)mesh.getSurfaceID(_triangle_3).values[i];
  dumper_surface.AddElemDataField(surf_id, 1, "surface_id");
  dumper_surface.SetPrefix("paraview/");
  dumper_surface.Init();
  dumper_surface.Dump();

  delete [] surf_id;
#endif //AKANTU_USE_IOHELPER

  finalize();
  return EXIT_SUCCESS;
}
