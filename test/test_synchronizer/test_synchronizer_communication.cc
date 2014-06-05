/**
 * @file   test_synchronizer_communication.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sun Sep 12 23:37:43 2010
 *
 * @brief  test to synchronize barycenters
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
#include "static_communicator.hh"
#include "distributed_synchronizer.hh"
#include "synchronizer_registry.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "mesh_partition_scotch.hh"
#include "communication_buffer.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "dumper_paraview.hh"
#endif //AKANTU_USE_IOHELPER

#include "test_data_accessor.hh"


using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);

  int dim = 2;
  ElementType type = _triangle_3;

  Mesh mesh(dim);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  bool wait=true;
  if(argc >1) {
    if(prank == 0)
    while(wait);
  }

  DistributedSynchronizer * communicator;
  if(prank == 0) {
    MeshIOMSH mesh_io;
    mesh_io.read("cube.msh", mesh);

    Mesh::type_iterator it = mesh.firstType();
    Mesh::type_iterator last_type = mesh.lastType();

    ElementTypeMapArray<Real> barycenters("", "", 0);
    mesh.initElementTypeMapArray(barycenters, dim, dim);

    GhostType ghost_type = _not_ghost;

    for(; it != last_type; ++it) {
      Array<Real> & mesh_data_array = *mesh.getDataPointer<Real>("barycenters", *it, ghost_type, dim);
      Array<Real>::iterator< Vector<Real> > mesh_data_array_it = mesh_data_array.begin(dim);

      UInt nb_element = mesh.getNbElement(*it, ghost_type);
      Array<Real> & barycenter = barycenters(*it, ghost_type);
      barycenter.resize(nb_element);
      Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(dim);

      for (UInt elem = 0; elem < nb_element; ++elem) {
        mesh.getBarycenter(elem, *it, bary_it->storage(), ghost_type);
        mesh.getBarycenter(elem, *it, mesh_data_array_it->storage(), ghost_type);
        ++bary_it;
        ++mesh_data_array_it;
      }
      debug::setDebugLevel(dblDump);
      std::cout << "Mesh Data barycenters (type "<< *it << ") :" << std::endl;
      std::cout << mesh_data_array;
      debug::setDebugLevel(dblInfo);
    }

    MeshPartition * partition = new MeshPartitionScotch(mesh, dim);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;

  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  // Checking if the barycenters of the partitioned elements match the ones in the partitioned MeshData
  ElementTypeMapArray<Real> barycenters("", "", 0);
  mesh.initElementTypeMapArray(barycenters, dim, dim);


  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator tit = mesh.firstType(dim, ghost_type);
    Mesh::type_iterator last_type = mesh.lastType(dim, ghost_type);

    for(; tit != last_type; ++tit) {
      Array<Real> & mesh_data_array = *mesh.getDataPointer<Real>("barycenters", *tit, ghost_type, dim);
      Array<Real>::iterator< Vector<Real> > mesh_data_array_it = mesh_data_array.begin(dim);

      UInt nb_element = mesh.getNbElement(*tit, ghost_type);
      Array<Real> & barycenter = barycenters(*tit, ghost_type);
      barycenter.resize(nb_element);
      Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(dim);
      UInt nb_component = barycenter.getNbComponent();
      for (UInt elem = 0; elem < nb_element; ++elem) {
	mesh.getBarycenter(elem, *tit, bary_it->storage(), ghost_type);
	for(UInt k(0); k < nb_component; ++k) {
	  AKANTU_DEBUG_ASSERT(bary_it->operator()(k) == mesh_data_array_it->operator()(k), "Barycenter doesn't match the value in MeshData. Calculated one is: " << *bary_it << " while Mesh Data has: " << *mesh_data_array_it);
	}
	++bary_it;
	++mesh_data_array_it;
      }
      debug::setDebugLevel(dblTest);
      std::cout << "Mesh Data barycenters (type "<< *tit << ") :" << std::endl;
      std::cout << mesh_data_array;
      debug::setDebugLevel(dblInfo);
    }
  }

  AKANTU_DEBUG_INFO("Creating TestAccessor");
  TestAccessor test_accessor(mesh, barycenters);
  SynchronizerRegistry synch_registry(test_accessor);
  synch_registry.registerSynchronizer(*communicator,_gst_test);

  AKANTU_DEBUG_INFO("Synchronizing tag");
  synch_registry.synchronize(_gst_test);


  // Checking the tags in MeshData (not a very good test because they're all identical,
  // but still...)
Array<UInt> & tags = mesh.getData<UInt>("tag_0", type);
  Array<UInt>::const_iterator< Vector<UInt> > tags_it = tags.begin(1);
  Array<UInt>::const_iterator< Vector<UInt> > tags_end = tags.end(1);
  AKANTU_DEBUG_ASSERT(mesh.getNbElement(type) == tags.getSize(),
                      "The number of tags does not match the number of elements on rank " << prank << ".");
  std::cout << " I am rank " << prank << " and here's my MeshData dump (it should contain " << mesh.getNbElement(type) << " elements and it has " << tags.getSize() << "!) :" << std::endl;
  std::cout << std::hex;
  debug::setDebugLevel(dblTest);
  for(; tags_it != tags_end; ++tags_it) {
    //AKANTU_DEBUG_ASSERT(*tags_it == 1, "The tag does not match the expected value on rank " << prank << " (got " << *tags_it << " instead.");
    std::cout << tags_it->operator()(0) << " ";
  }
  debug::setDebugLevel(dblInfo);
  std::cout << std::endl;

#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper("test-scotch-partition");
  dumper.registerMesh(mesh, _all_dimensions, _not_ghost);
  dumper.registerField("partitions", new DumperIOHelper::ElementPartitionField<>(mesh, _all_dimensions, _not_ghost));
  dumper.dump();

  DumperParaview dumper_ghost("test-scotch-partition-ghost");
  dumper_ghost.registerMesh(mesh, _all_dimensions, _ghost);
  dumper_ghost.registerField("partitions", new DumperIOHelper::ElementPartitionField<>(mesh, _all_dimensions, _ghost));
  dumper_ghost.dump();

#endif //AKANTU_USE_IOHELPER
  delete communicator;
  finalize();

  return EXIT_SUCCESS;
}
