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
#include "synchronizer.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
GhostSynchronizer::GhostSynchronizer() :
  nb_ghost_synchronization_tags(0) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
GhostSynchronizer::~GhostSynchronizer() {
  AKANTU_DEBUG_IN();

  for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
       it != synchronizers.end();
       ++it) {
    delete (*it);
  }
  synchronizers.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::synchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(registered_synchronization.find(tag) != registered_synchronization.end(),
		      "Tag " << tag << " is not registered.");

  AKANTU_DEBUG_INFO("Synchronizing the tag : "
		    << registered_synchronization.find(tag)->second
		    << " (" << tag <<")");
  for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
       it != synchronizers.end();
       ++it) {
    (*it)->synchronize(tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::asynchronousSynchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(registered_synchronization.find(tag) != registered_synchronization.end(),
		      "Tag " << tag << " is not registered.");

  for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
       it != synchronizers.end();
       ++it) {
    (*it)->asynchronousSynchronize(tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::waitEndSynchronize(GhostSynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(registered_synchronization.find(tag) != registered_synchronization.end(),
		      "Tag " << tag << " is not registered.");

  for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
       it != synchronizers.end();
       ++it) {
    (*it)->waitEndSynchronize(tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::registerTag(GhostSynchronizationTag tag,
				    const std::string & name) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(registered_synchronization.find(tag) == registered_synchronization.end(),
		      "Tag " << tag << " (" << name << ") already registered.");

  registered_synchronization[tag] = name;
  for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
       it != synchronizers.end();
       ++it) {
    (*it)->registerTag(tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GhostSynchronizer::registerSynchronizer(Synchronizer & synchronizer) {
  AKANTU_DEBUG_IN();

  synchronizers.push_back(&synchronizer);

  synchronizer.registerGhostSynchronizer(*this);

  std::map<GhostSynchronizationTag, std::string>::iterator it;
  for (it = registered_synchronization.begin();
       it != registered_synchronization.end();
       ++it) {
    synchronizer.registerTag(it->first);
  }


  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
