/**
 * @file   aka_common.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Jun 11 16:56:43 2010
 *
 * @brief Initialization of global variables
 *
 * @section LICENSE
 *
 * <insert license here>
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
