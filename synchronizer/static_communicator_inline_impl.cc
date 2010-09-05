/**
 * @file   static_communicator_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Sep  6 00:16:19 2010
 *
 * @brief  implementation of inline functions
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::freeCommunicationRequest(CommunicationRequest * request) {
  delete request;
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::freeCommunicationRequest(std::vector<CommunicationRequest *> & requests) {
  std::vector<CommunicationRequest *>::iterator it;
  for(it = requests.begin(); it != requests.end(); ++it) {
    delete (*it);
  }
}
