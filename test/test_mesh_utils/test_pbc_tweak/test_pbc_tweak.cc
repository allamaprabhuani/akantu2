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
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);
  
  MeshUtils::tweakConnectivityForPBC(mesh,1,0,0);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  model->initVectors();
  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 2*nb_nodes*sizeof(Real));
  memset(model->getVelocity().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().values, 0, 2*nb_nodes*sizeof(Real));

  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

 #ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-pbc-tweak");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_tetrahedron_4).values,
			 TETRA1, mesh.getNbElement(_tetrahedron_4), C_MODE);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}
