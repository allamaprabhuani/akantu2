/**
 * @file   aka_common.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 16:56:43 2010
 *
 * @brief Initialization of global variables
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
#include "aka_common.hh"
#include "aka_static_memory.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void initialize(int * argc, char *** argv) {
  AKANTU_DEBUG_IN();

  StaticMemory::getStaticMemory();
  StaticCommunicator * comm = StaticCommunicator::getStaticCommunicator(argc, argv);
  debug::setParallelContext(comm->whoAmI(), comm->getNbProc());
  debug::initSignalHandler();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void finalize() {
  AKANTU_DEBUG_IN();

  if(StaticMemory::isInstantiated()) delete StaticMemory::getStaticMemory();
  if(StaticCommunicator::isInstantiated()) {
    StaticCommunicator *comm = StaticCommunicator::getStaticCommunicator();
    comm->barrier();
    delete comm;
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
