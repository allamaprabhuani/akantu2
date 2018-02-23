/**
 * @file   static_communicator_dummy.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 13 2016
 *
 * @brief  Dummy communicator to make everything work im sequential
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <type_traits>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

Communicator::Communicator(int & /*argc*/, char **& /*argv*/,
                           const private_member & /*unused*/) {}

template <typename T>
void Communicator::sendImpl(const T *, Int, Int, Int,
                            const CommunicationMode &) const {}
template <typename T>
void Communicator::receiveImpl(T *, Int, Int, Int) const {}

template <typename T>
CommunicationRequest
Communicator::asyncSendImpl(const T *, Int, Int, Int,
                            const CommunicationMode &) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
CommunicationRequest Communicator::asyncReceiveImpl(T *, Int, Int, Int) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::probe(Int, Int, CommunicationStatus &) const {}
template <typename T>
bool Communicator::asyncProbe(Int, Int, CommunicationStatus &) const {
  return true;
}

bool Communicator::test(CommunicationRequest &) const { return true; }
bool Communicator::testAll(std::vector<CommunicationRequest> &) const {
  return true;
}
void Communicator::wait(CommunicationRequest &) const {}
void Communicator::waitAll(std::vector<CommunicationRequest> &) const {}
UInt Communicator::waitAny(std::vector<CommunicationRequest> &) const {
  return UInt(-1);
}

void Communicator::barrier() const {}
CommunicationRequest Communicator::asyncBarrier() const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::reduceImpl(T *, int, const SynchronizerOperation &,
                              int) const {}

template <typename T>
void Communicator::allReduceImpl(T *, int,
                                 const SynchronizerOperation &) const {}

template <typename T> inline void Communicator::allGatherImpl(T *, int) const {}
template <typename T>
inline void Communicator::allGatherVImpl(T *, int *) const {}

template <typename T>
inline void Communicator::gatherImpl(T *, int, int) const {}
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, T * gathered,
                              int) const {
  static_assert(std::is_trivially_copyable<T>{},
                "Cannot send this type of data");
  std::memcpy(gathered, values, nb_values);
}

template <typename T>
inline void Communicator::gatherVImpl(T *, int *, int) const {}
template <typename T>
inline void Communicator::broadcastImpl(T *, int, int) const {}

int Communicator::getMaxTag() const { return std::numeric_limits<int>::max(); }
int Communicator::getMinTag() const { return 0; }

} // namespace akantu
