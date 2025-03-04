/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <cstring>
#include <type_traits>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

// NOLINTBEGIN(misc-definitions-in-headers)

Communicator::Communicator(int & /*argc*/, char **& /*argv*/,
                           const private_member & /*unused*/) {}

Communicator::Communicator(const private_member & /*unused*/) {}

template <typename T>
void Communicator::sendImpl(const T * /*unused*/, Int /*unused*/,
                            Int /*unused*/, Int /*unused*/,
                            const CommunicationMode & /*unused*/) const {}
template <typename T>
void Communicator::receiveImpl(T * /*unused*/, Int /*unused*/, Int /*unused*/,
                               Int /*unused*/) const {}

template <typename T>
CommunicationRequest
Communicator::asyncSendImpl(const T * /*unused*/, Int /*unused*/,
                            Int /*unused*/, Int /*unused*/,
                            const CommunicationMode & /*unused*/) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
CommunicationRequest
Communicator::asyncReceiveImpl(T * /*unused*/, Int /*unused*/, Int /*unused*/,
                               Int /*unused*/) const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::probe(Int /*unused*/, Int /*unused*/,
                         CommunicationStatus & /*unused*/) const {}
template <typename T>
bool Communicator::asyncProbe(Int /*unused*/, Int /*unused*/,
                              CommunicationStatus & /*unused*/) const {
  return true;
}

bool Communicator::test(CommunicationRequest & /*unused*/) { return true; }
bool Communicator::testAll(std::vector<CommunicationRequest> & /*unused*/) {
  return true;
}
void Communicator::wait(CommunicationRequest & /*unused*/) {}
void Communicator::waitAll(std::vector<CommunicationRequest> & /*unused*/) {}
Int Communicator::waitAny(std::vector<CommunicationRequest> & /*unused*/) {
  return -1;
}

void Communicator::barrier() const {}
CommunicationRequest Communicator::asyncBarrier() const {
  return std::shared_ptr<InternalCommunicationRequest>(
      new InternalCommunicationRequest(0, 0));
}

template <typename T>
void Communicator::reduceImpl(T * /*unused*/, int /*unused*/,
                              SynchronizerOperation /*unused*/,
                              int /*unused*/) const {}

template <typename T>
void Communicator::allReduceImpl(T * /*unused*/, int /*unused*/,
                                 SynchronizerOperation /*unused*/) const {}

template <typename T>
void Communicator::scanImpl(T * values, T * result, int n,
                            SynchronizerOperation /*unused*/) const {
  if (values == result) {
    return;
  }

  std::copy_n(values, n, result);
}

template <typename T>
void Communicator::exclusiveScanImpl(T * /*values*/, T * result, int n,
                                     SynchronizerOperation /*unused*/) const {
  std::fill_n(result, n, T());
}

template <typename T>
void Communicator::allGatherImpl(T * /*unused*/, int /*unused*/) const {}
template <typename T>
void Communicator::allGatherVImpl(T * /*unused*/,
                                  const int * /*unused*/) const {}

template <typename T>
void Communicator::gatherImpl(T * /*unused*/, int /*unused*/,
                              int /*unused*/) const {}
template <typename T>
void Communicator::gatherImpl(T * values, int nb_values, T * gathered,
                              int /*unused*/) const {
  static_assert(std::is_trivially_copyable<T>{},
                "Cannot send this type of data");
  std::memcpy(gathered, values, nb_values);
}

template <typename T>
void Communicator::gatherVImpl(T * /*unused*/, int * /*unused*/,
                               int /*unused*/) const {}
template <typename T>
void Communicator::broadcastImpl(T * /*unused*/, int /*unused*/,
                                 int /*unused*/) const {}

int Communicator::getMaxTag() const { return std::numeric_limits<int>::max(); }
int Communicator::getMinTag() const { return 0; }

Int Communicator::getNbProc() const { return 1; }
Int Communicator::whoAmI() const { return 0; }

// NOLINTEND(misc-definitions-in-headers)

} // namespace akantu
