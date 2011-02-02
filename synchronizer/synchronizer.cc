/**
 * @file   synchronizer.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 26 16:31:48 2010
 *
 * @brief  implementation of the common part
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
