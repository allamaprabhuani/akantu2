/**
 * @file   distributed_synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Thu Aug 19 15:28:35 2010
 *
 * @brief  wrapper to the static communicator
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

#ifndef __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__
#define __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "static_communicator.hh"
#include "synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class CommunicationBuffer;

class DistributedSynchronizer : public Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DistributedSynchronizer(SynchronizerID id = "distributedSynchronizer", 
			  MemoryID memory_id = 0);
  virtual ~DistributedSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a mesh and a partition and create the local mesh and the associated distributedSynchronizer
  static DistributedSynchronizer * 
  createDistributedSynchronizerMesh(Mesh & mesh,
				    const MeshPartition * partition,
				    UInt root = 0,
				    SynchronizerID id = "distributedSynchronizer",
				    MemoryID memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

protected:
  /// fill the nodes type vector
  void fillNodesType(Mesh & mesh);


  /// fill the communications array of a distributedSynchronizer based on a partition array
  void fillCommunicationScheme(UInt * partition,
			       UInt nb_local_element,
			       UInt nb_ghost_element,
			       ElementType type);

  /// compute buffer size for a given tag and data accessor 
  void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the static memory instance
  StaticCommunicator * static_communicator;

  /// size of data to send to each processor by communication tag
  std::map< SynchronizationTag, Vector<UInt> > size_to_send;

  /// size of data to receive form each processor by communication tag
  std::map< SynchronizationTag, Vector<UInt> > size_to_receive;

  CommunicationBuffer * send_buffer;
  CommunicationBuffer * recv_buffer;

  /// send requests
  std::vector<CommunicationRequest *> send_requests;
  /// receive requests
  std::vector<CommunicationRequest *> recv_requests;

  /// list of real element to send ordered by type then by receiver processors
  ByElementTypeUInt element_to_send_offset;
  ByElementTypeUInt element_to_send;

  std::vector<Element> * send_element;
  std::vector<Element> * recv_element;

  /// list of ghost element to receive ordered by type then by sender processors
  ByElementTypeUInt element_to_receive_offset;
  ByElementTypeUInt element_to_receive;

  UInt nb_proc;
  UInt rank;

  /// global node ids
  Vector<UInt> * nodes_global_ids;

  /// node type,  -3 pure ghost, -2  master for the  node, -1 normal node,  i in
  /// [0-N] slave node and master is proc i
  Vector<Int> * nodes_type;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "distributedSynchronizer_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__ */
