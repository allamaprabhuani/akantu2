/**
 * @file   synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 26 16:31:48 2010
 *
 * @brief  implementation of the common part
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "synchronizer.hh"
#include "ghost_synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Synchronizer::Synchronizer(SynchronizerID id, MemoryID memory_id) :
  Memory(memory_id), id(id) {

}

/* -------------------------------------------------------------------------- */
void Synchronizer::registerGhostSynchronizer(const GhostSynchronizer & ghost_synchronizer) {
  AKANTU_DEBUG_IN();
  this->ghost_synchronizer = &ghost_synchronizer;
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
