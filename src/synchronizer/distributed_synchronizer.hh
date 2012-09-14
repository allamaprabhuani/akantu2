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

class DistributedSynchronizer : public Synchronizer, public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  DistributedSynchronizer(Mesh & mesh,
			  SynchronizerID id = "distributed_synchronizer",
			  MemoryID memory_id = 0);

public:
  virtual ~DistributedSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a  mesh and a partition and  create the local mesh  and the associated
  /// DistributedSynchronizer
  static DistributedSynchronizer *
  createDistributedSynchronizerMesh(Mesh & mesh,
				    const MeshPartition * partition,
				    UInt root = 0,
				    SynchronizerID id = "distributed_synchronizer",
				    MemoryID memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// mesh event handler onRemovedElement
  virtual void onElementsRemoved(const Vector<Element> & element_list,
				 const ByElementTypeUInt & new_numbering);

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
protected:
  /// reference to the underlying mesh
  Mesh & mesh;

  /// the static memory instance
  StaticCommunicator * static_communicator;

  class Communication {
  public:
    void resize(UInt size) {
      send_buffer.resize(size);
      recv_buffer.resize(size);
      size_to_send   .resize(size);
      size_to_receive.resize(size);
    }

  public:
    /// size of data to send to each processor
    std::vector<UInt> size_to_send;
    /// size of data to recv to each processor
    std::vector<UInt> size_to_receive;
    std::vector<CommunicationBuffer> send_buffer;
    std::vector<CommunicationBuffer> recv_buffer;

    std::vector<CommunicationRequest *> send_requests;
    std::vector<CommunicationRequest *> recv_requests;
  };

  std::map<SynchronizationTag, Communication> communications;

  /// list of element to sent to proc p
  Vector<Element> * send_element;
  /// list of element to receive from proc p
  Vector<Element> * recv_element;

  UInt nb_proc;
  UInt rank;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "distributedSynchronizer_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_DISTRIBUTED_SYNCHRONIZER_HH__ */
