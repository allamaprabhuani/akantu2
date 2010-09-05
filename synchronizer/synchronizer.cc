/**
 * @file   synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 26 16:31:48 2010
 *
 * @brief  implementation of the common part
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "synchronizer.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

Synchronizer::Synchronizer(SynchronizerID id, MemoryID memory_id) :
  Memory(memory_id), id(id) {

}

__END_AKANTU__
