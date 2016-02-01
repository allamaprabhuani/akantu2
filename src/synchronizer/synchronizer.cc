/**
 * @file   synchronizer.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 01 2010
 * @date last modification: Wed Nov 13 2013
 *
 * @brief  implementation of the common part
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Synchronizer::Synchronizer(SynchronizerID id, MemoryID memory_id,
			   StaticCommunicator & comm) :
  Memory(id, memory_id),
  static_communicator(&comm) {

}

/* -------------------------------------------------------------------------- */
void Synchronizer::synchronize(DataAccessor & data_accessor,
			       SynchronizationTag tag) {
  AKANTU_DEBUG_IN();
  asynchronousSynchronize(data_accessor,tag);
  waitEndSynchronize(data_accessor,tag);
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */


__END_AKANTU__
