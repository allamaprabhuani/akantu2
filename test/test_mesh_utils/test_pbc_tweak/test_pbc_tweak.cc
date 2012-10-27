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
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  initialize(argc, argv);
  debug::setDebugLevel(akantu::dblInfo);

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);

  SolidMechanicsModel model(mesh);
  /* -------------------------------------------------------------------------- */
  model.initFull("material.dat");
  /* -------------------------------------------------------------------------- */
  model.setPBC(1,1,1);
  model.initPBC();
  model.assembleMassLumped();
  /* -------------------------------------------------------------------------- */

  model.setBaseName("test-pbc-tweak");
  model.addDumpField("mass");
  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
