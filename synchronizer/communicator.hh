/**
 * @file   communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:28:35 2010
 *
 * @brief  wrapper to the static communicator
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATOR_HH__
#define __AKANTU_COMMUNICATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "static_communicator.hh"
#include "synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Communicator : public Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Communicator(SynchronizerID id = "communicator", MemoryID memory_id = 0);
  virtual ~Communicator();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get a mesh and a partition and create the local mesh and the associated communicator
  static Communicator * createCommunicatorDistributeMesh(Mesh & mesh,
							 const MeshPartition * partition,
							 UInt root = 0,
							 SynchronizerID id = "communicator",
							 MemoryID memory_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Inherited from Synchronizer                                              */
  /* ------------------------------------------------------------------------ */
  /// synchronize ghosts
  void synchronize(GhostSynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(GhostSynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(GhostSynchronizationTag tag);

  /// do a all reduce operation
  void allReduce(Real * values, UInt nb_values, const SynchronizerOperation & op);

  /* ------------------------------------------------------------------------ */
  /// register a new communication
  void registerTag(GhostSynchronizationTag tag);

protected:
  /// fill the communications array of a communicator based on a partition array
  void fillCommunicationScheme(UInt * partition,
			       UInt nb_local_element,
			       UInt nb_ghost_element,
			       ElementType type);

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
  std::map< GhostSynchronizationTag, Vector<UInt> > size_to_send;

  /// size of data to receive form each processor by communication tag
  std::map< GhostSynchronizationTag, Vector<UInt> > size_to_receive;

  Vector<Real> * send_buffer;
  Vector<Real> * recv_buffer;

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
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "communicator_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_COMMUNICATOR_HH__ */
