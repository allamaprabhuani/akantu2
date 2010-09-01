/**
 * @file   ghost_synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Aug 24 18:59:17 2010
 *
 * @brief  implementation of the ghost synchronizer
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "ghost_synchronizer.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
GhostSynchronizer::GhostSynchronizer() :
  nb_ghost_synchronization_tags(0) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::registerTag(GhostSynchronizationTag tag,
				    const std::string & name) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(registered_synchronization.find(tag) == registered_synchronization.end(),
		      "Tag " << tag << " (" << name << ") already registered.");

  registered_synchronization[tag] = name;

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
