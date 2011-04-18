/**
 * @file   static_communicator_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Sep  6 00:16:19 2010
 *
 * @brief  implementation of inline functions
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

/* -------------------------------------------------------------------------- */
#define AKANTU_BOOST_REAL_COMMUNICATOR_CALL(r, call, comm_type)		\
  case BOOST_PP_LIST_AT(comm_type, 0): {				\
    BOOST_PP_LIST_AT(comm_type, 1) * comm =				\
      dynamic_cast<BOOST_PP_LIST_AT(comm_type, 1) *>(real_static_communicator); \
    BOOST_PP_IF(BOOST_PP_LIST_AT(call, 0),				\
		return comm->BOOST_PP_LIST_AT(call, 1),			\
		comm->BOOST_PP_LIST_AT(call, 1));			\
    break;								\
  }

#define AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(call, ret)		\
  do {									\
    switch(real_type)							\
      {									\
	BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_REAL_COMMUNICATOR_CALL,	\
			      (ret, (call, BOST_PP_NIL)),		\
			      AKANTU_COMMUNICATOR_LIST_ALL)		\
      default:								\
	{								\
	  AKANTU_DEBUG_ERROR("Wrong communicator : " << real_type);	\
	}								\
      }									\
  } while(0)


/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicator::send(T * buffer, Int size, Int receiver, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(send(buffer, size, receiver, tag), 0);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicator::receive(T * buffer, Int size, Int sender, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(receive(buffer, size, sender, tag), 0);
}


/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationRequest * StaticCommunicator::asyncSend(T * buffer, Int size,
							    Int receiver, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(asyncSend(buffer, size, receiver, tag), 1);
  return NULL;
}

/* -------------------------------------------------------------------------- */
template<typename T> inline CommunicationRequest * StaticCommunicator::asyncReceive(T * buffer, Int size,
										    Int sender, Int tag) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(asyncReceive(buffer, size, sender, tag), 1);
  return NULL;
}

/* -------------------------------------------------------------------------- */
template<typename T> inline void StaticCommunicator::allReduce(T * values, Int nb_values,
							       const SynchronizerOperation & op) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(allReduce(values, nb_values, op), 0);
}

/* -------------------------------------------------------------------------- */
template<typename T> inline void StaticCommunicator::gather(T * values, Int nb_values, Int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(gather(values, nb_values, root), 0);
}

/* -------------------------------------------------------------------------- */
template<typename T> inline void StaticCommunicator::gatherv(T * values, Int * nb_values, Int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(gatherv(values, nb_values, root), 0);
}

/* -------------------------------------------------------------------------- */
template<typename T> inline void StaticCommunicator::broadcast(T * values, Int nb_values, Int root) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(broadcast(values, nb_values, root), 0);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::barrier() {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(barrier(), 0);
}

/* -------------------------------------------------------------------------- */
inline bool StaticCommunicator::testRequest(CommunicationRequest * request) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(testRequest(request), 1);
  return false;
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::wait(CommunicationRequest * request) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(wait(request), 0);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicator::waitAll(std::vector<CommunicationRequest *> & requests) {
  AKANTU_BOOST_REAL_COMMUNICATOR_SELECT_CALL(waitAll(requests), 0);
}
