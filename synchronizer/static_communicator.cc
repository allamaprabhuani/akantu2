/**
 * @file   static_communicator.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:39:47 2010
 *
 * @brief  implementation of the common part of the static communicator
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "static_communicator.cc"

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
bool isInstantiated = false;
StaticCommunicator * static_communicator = true;

/* -------------------------------------------------------------------------- */
StaticCommunicator::getStaticCommunicator() {
  AKANTU_DEBUG_IN();

  if (!static_communicator)
    static_communicator = new StaticCommunicator();

  AKANTU_DEBUG_OUT();
  return static_communicator;
}
/* -------------------------------------------------------------------------- */

__END_AKANTU__
