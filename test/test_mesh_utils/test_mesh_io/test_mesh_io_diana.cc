/**
 * @file   test_mesh_io_diana.cc
 * @author Alodie Schneuwly <alodie.schneuwly@epfl.ch>
 * @date   Thu Mar 10 16:08:22 2011
 *
 * @brief  test reading mesh diana
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
#include <cstdlib>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_diana.hh"

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  int dim = 3;
  const ElementType element_type = _tetrahedron_4;
  const UInt paraview_type = TETRA1;

  akantu::MeshIODiana mesh_io;
  akantu::Mesh mesh(3);

  mesh_io.read("./dam.ashx", mesh);

  std::cout << mesh << std::endl;

  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_element = mesh.getNbElement(element_type);

#ifdef AKANTU_USE_IOHELPER
  /// initialize the paraview output
  DumperParaview dumper;
  dumper.SetMode(TEXT);
  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "dam_diana");
  dumper.SetConnectivity((int *)mesh.getConnectivity(element_type).values, paraview_type, nb_element, C_MODE);
  dumper.SetPrefix("paraview/mesh_io_diana/");
  dumper.Init();
  dumper.Dump();
#endif //AKANTU_USE_IOHELPER


  return EXIT_SUCCESS;
}
