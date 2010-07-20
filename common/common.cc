/**
 * @file   common.cc
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
#include "common.hh"
#include "static_memory.hh"

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
  delete StaticMemory::getStaticMemory();
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
