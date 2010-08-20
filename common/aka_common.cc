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
#ifdef AKANTU_USE_MPI
#  include "static_communicator.hh"
#endif
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void initialize() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void finalize() {
  AKANTU_DEBUG_IN();

  if(StaticMemory::isInstantiated()) delete StaticMemory::getStaticMemory();

#ifdef AKANTU_USE_MPI
  if(StaticCommunicator::isInstantiated()) delete StaticCommunicator::getStaticMemory();
#endif

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
