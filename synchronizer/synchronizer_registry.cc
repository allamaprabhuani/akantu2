/**
 * @file   synchronizer_registry.cc
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Wed Jun 15 15:56:07 2011
 *
 * @brief  Registry of synchronizers
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
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */
SynchronizerRegistry::SynchronizerRegistry(DataAccessor & da) :
  // nb_synchronization_tags(0),
  data_accessor(da) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SynchronizerRegistry::~SynchronizerRegistry() {
  AKANTU_DEBUG_IN();

  // for (std::list<Synchronizer *>::iterator it = synchronizers.begin();
  //      it != synchronizers.end();
  //      ++it) {
  //   delete (*it);
  // }
  // synchronizers.clear();

  synchronizers.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::synchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  std::pair<Tag2Sync::iterator,Tag2Sync::iterator> range =
    synchronizers.equal_range(tag);

  for (Tag2Sync::iterator it = range.first; it != range.second;++it) {
    (*it).second->synchronize(data_accessor,tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::asynchronousSynchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  std::pair<Tag2Sync::iterator,Tag2Sync::iterator> range =
    synchronizers.equal_range(tag);

  for (Tag2Sync::iterator it = range.first; it != range.second;++it) {
    (*it).second->asynchronousSynchronize(data_accessor,tag);
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::waitEndSynchronize(SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  std::pair<Tag2Sync::iterator,Tag2Sync::iterator> range =
    synchronizers.equal_range(tag);

  for (Tag2Sync::iterator it = range.first; it != range.second;++it) {
    (*it).second->waitEndSynchronize(data_accessor,tag);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SynchronizerRegistry::registerSynchronizer(Synchronizer & synchronizer,
						SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  synchronizers.
    insert(std::pair<SynchronizationTag,Synchronizer *>(tag,&synchronizer));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
__END_AKANTU__
