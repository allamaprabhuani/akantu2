/**
 * @file   synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 23 13:48:37 2010
 *
 * @brief  interface for communicator and pbc synchronizers
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_HH__
#define __AKANTU_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {
  class GhostSynchronizer;
}

__BEGIN_AKANTU__

class Synchronizer : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Synchronizer(SynchronizerID id = "synchronizer", MemoryID memory_id = 0);

  virtual ~Synchronizer() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// synchronize ghosts
  virtual void synchronize(GhostSynchronizationTag tag) = 0;

  /// asynchronous synchronization of ghosts
  virtual void asynchronousSynchronize(GhostSynchronizationTag tag) = 0;

  /// wait end of asynchronous synchronization of ghosts
  virtual void waitEndSynchronize(GhostSynchronizationTag tag) = 0;

  virtual void allReduce(Real * values, UInt nb_values, const SynchronizerOperation & op) = 0;

  /* ------------------------------------------------------------------------ */
  /// register a new communication
  virtual void registerTag(GhostSynchronizationTag tag) = 0;

  void registerGhostSynchronizer(const GhostSynchronizer & ghost_synchronizer);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// id of the synchronizer
  SynchronizerID id;

  /// the associated ghost_synchonizer
  const GhostSynchronizer * ghost_synchronizer;
};


__END_AKANTU__

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
