/**
 * @file   communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:28:35 2010
 *
 * @brief  wrapper to the static communicator
 *
 * @section LICENSE
 *
 * <insert license here>
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

class Communicator : Synchronizer {
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

  /// fill the communications array of a communicator based on a partition array
  void fillCommunicationScheme(UInt * partition,
			       UInt nb_local_element,
			       UInt nb_ghost_element,
			       UInt nb_element_to_send,
			       ElementType type);

protected:

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

  Vector<Real> send_buffer;
  Vector<Real> receive_buffer;

  /// list of real element to send ordered by type then by receiver processors
  ByElementTypeUInt element_to_send_offset;
  ByElementTypeUInt element_to_send;

  /// list of ghost element to receive ordered by type then by sender processors
  ByElementTypeUInt element_to_receive_offset;
  ByElementTypeUInt element_to_receive;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "communicator_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_COMMUNICATOR_HH__ */
