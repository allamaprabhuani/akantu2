/**
 * @file   test_mesh_partitionate_scotch.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Sep 01 17:57:12 2010
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
#include "mesh_partition_scotch.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;
/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);
  debug::setDebugLevel(akantu::dblDump);

  int dim = 2;

  akantu::Mesh mesh(dim);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);

  akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh, dim);
  partition->partitionate(8);

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper("test-scotch-partition");
  DumperIOHelper::Field * field = new DumperIOHelper::ElementalField<UInt>(partition->getPartitions(),
									   dim);
  dumper.registerMesh(mesh, dim);
  dumper.registerField("partitions", field);
  dumper.dump();
#endif //AKANTU_USE_IOHELPER

  partition->reorder();
  mesh_io.write("triangle_reorder.msh", mesh);

  delete partition;

  finalize();

  return EXIT_SUCCESS;
}
