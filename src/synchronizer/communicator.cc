/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#if defined(AKANTU_USE_MPI)
#include "mpi_communicator_data.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

Int InternalCommunicationRequest::counter = 0;

/* -------------------------------------------------------------------------- */
InternalCommunicationRequest::InternalCommunicationRequest(Idx source, Idx dest)
    : source(source), destination(dest) {
  this->id = counter++;
}

/* -------------------------------------------------------------------------- */
InternalCommunicationRequest::~InternalCommunicationRequest() = default;

/* -------------------------------------------------------------------------- */
void InternalCommunicationRequest::printself(std::ostream & stream,
                                             int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "CommunicationRequest [" << std::endl;
  stream << space << " + id          : " << id << std::endl;
  stream << space << " + source      : " << source << std::endl;
  stream << space << " + destination : " << destination << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
Communicator::~Communicator() {
  auto * event = new FinalizeCommunicatorEvent(*this);
  this->sendEvent(*event);
  delete event;
}

/* -------------------------------------------------------------------------- */
Communicator & Communicator::getStaticCommunicator() {
  return getWorldCommunicator();
}

/* -------------------------------------------------------------------------- */
Communicator & Communicator::getStaticCommunicator(int & /*argc*/,
                                                   char **& /*argv*/) {
  return getWorldCommunicator();
}

} // namespace akantu

#ifdef AKANTU_USE_MPI
#include "communicator_mpi_inline_impl.hh"
#else
#include "communicator_dummy_inline_impl.hh"
#endif

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Template instantiation                                                     */
/* -------------------------------------------------------------------------- */
#define AKANTU_COMM_INSTANTIATE(T)                                             \
  template void Communicator::probe<T>(Int sender, Int tag,                    \
                                       CommunicationStatus & status) const;    \
  template bool Communicator::asyncProbe<T>(                                   \
      Int sender, Int tag, CommunicationStatus & status) const;                \
  template void Communicator::sendImpl<T>(                                     \
      const T * buffer /*NOLINT*/, Int size, Int receiver, Int tag,            \
      const CommunicationMode & mode) const;                                   \
  template void Communicator::receiveImpl<T>(T * buffer /*NOLINT*/, Int size,  \
                                             Int sender, Int tag) const;       \
  template CommunicationRequest Communicator::asyncSendImpl<T>(                \
      const T * buffer /*NOLINT*/, Int size, Int receiver, Int tag,            \
      const CommunicationMode & mode) const;                                   \
  template CommunicationRequest Communicator::asyncReceiveImpl<T>(             \
      T * buffer /* NOLINT */, Int size, Int sender, Int tag) const;           \
  template void Communicator::allGatherImpl<T>(T * values /*NOLINT*/,          \
                                               int nb_values) const;           \
  template void Communicator::allGatherVImpl<T>(T * values /*NOLINT*/,         \
                                                const int * nb_values) const;  \
  template void Communicator::gatherImpl<T>(T * values /*NOLINT*/,             \
                                            int nb_values, int root) const;    \
  template void Communicator::gatherImpl<T>(                                   \
      T * values /*NOLINT*/, int nb_values, T * gathered /*NOLINT*/,           \
      int nb_gathered) const;                                                  \
  template void Communicator::gatherVImpl<T>(T * values /*NOLINT*/,            \
                                             int * nb_values, int root) const; \
  template void Communicator::broadcastImpl<T>(T * values /*NOLINT*/,          \
                                               int nb_values, int root) const; \
  template void Communicator::allReduceImpl<T>(                                \
      T * values /*NOLINT*/, int nb_values, SynchronizerOperation op) const;   \
  template void Communicator::scanImpl<T>(T * values /*NOLINT*/,               \
                                          T * /*NOLINT*/, int nb_values,       \
                                          SynchronizerOperation op) const;     \
  template void Communicator::exclusiveScanImpl<T>(                            \
      T * values /*NOLINT*/, T * /*NOLINT*/, int nb_values,                    \
      SynchronizerOperation op) const

#define MIN_MAX_REAL SCMinMaxLoc<Real, int>

#if !defined(DOXYGEN)

AKANTU_COMM_INSTANTIATE(bool);
AKANTU_COMM_INSTANTIATE(float);
AKANTU_COMM_INSTANTIATE(double);
AKANTU_COMM_INSTANTIATE(unsigned int);
AKANTU_COMM_INSTANTIATE(long long);
AKANTU_COMM_INSTANTIATE(unsigned long long);
AKANTU_COMM_INSTANTIATE(int);
AKANTU_COMM_INSTANTIATE(char);
AKANTU_COMM_INSTANTIATE(NodeFlag);
AKANTU_COMM_INSTANTIATE(MIN_MAX_REAL);
AKANTU_COMM_INSTANTIATE(std::size_t);

#if AKANTU_INTEGER_SIZE > 4
AKANTU_COMM_INSTANTIATE(int);
#endif

#endif

} // namespace akantu
