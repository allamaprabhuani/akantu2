/**
 * @file   communications_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Tue Sep  6 17:14:06 2016
 *
 * @brief  Implementation of Communications
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
#include "communications.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATIONS_TMPL_HH__
#define __AKANTU_COMMUNICATIONS_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications<Entity>::iterator {
  typedef
      typename std::map<UInt, Communication>::iterator communication_iterator;

public:
  iterator() : communications(NULL) {}
  iterator(scheme_iterator scheme_it, communication_iterator comm_it,
           Communications<Entity> & communications, const SynchronizationTag & tag)
      : scheme_it(scheme_it), comm_it(comm_it), communications(&communications),
        tag(tag) {}

  iterator & operator=(const iterator & other) {
    if (this != &other) {
      this->scheme_it = other.scheme_it;
      this->comm_it = other.comm_it;
      this->communications = other.communications;
      this->tag = other.tag;
    }
    return *this;
  }

  iterator & operator++() {
    ++scheme_it;
    ++comm_it;
    return *this;
  }

  CommunicationDescriptor<Entity> operator*() {
    AKANTU_DEBUG_ASSERT(
        scheme_it->first == comm_it->first,
        "The to iterators are not in phase, something wrong"
            << " happened, time to take out your favorite debugger");
    return CommunicationDescriptor<Entity>(comm_it->second, scheme_it->second,
                                           *communications, tag,
                                           scheme_it->first);
  }

  bool operator==(const iterator & other) const {
    return scheme_it == other.scheme_it && comm_it == other.comm_it;
  }

  bool operator!=(const iterator & other) const {
    return scheme_it != other.scheme_it || comm_it != other.comm_it;
  }

private:
  scheme_iterator scheme_it;
  communication_iterator comm_it;
  Communications<Entity> * communications;
  SynchronizationTag tag;
};

/* -------------------------------------------------------------------------- */
template <class Entity> class Communications<Entity>::tag_iterator {
  typedef std::map<SynchronizationTag, UInt>::const_iterator internal_iterator;

public:
  tag_iterator(const internal_iterator & it) : it(it) {}
  tag_iterator & operator++() {
    ++it;
    return *this;
  }
  SynchronizationTag operator*() { return it->first; }
  bool operator==(const tag_iterator & other) const { return it == other.it; }
  bool operator!=(const tag_iterator & other) const { return it != other.it; }

private:
  internal_iterator it;
};

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::CommunicationPerProcs &
Communications<Entity>::getCommunications(const SynchronizationTag & tag,
                                          const CommunicationSendRecv & sr) {
  CommunicationsPerTags::iterator comm_it = communications[sr].find(tag);
  if (comm_it == communications[sr].end())
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "No known communications for the tag: " << tag);
  return comm_it->second;
}

/* ---------------------------------------------------------------------- */
template <class Entity>
UInt Communications<Entity>::getPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) const {
  const std::map<SynchronizationTag, UInt> & pending =
      pending_communications[sr];
  std::map<SynchronizationTag, UInt>::const_iterator it = pending.find(tag);
  AKANTU_DEBUG_ASSERT(it != pending.end(),
                      "Pending communications are not initialized!");
  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
bool Communications<Entity>::hasPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) const {
  return this->hasCommunication(tag) && (this->getPending(tag, sr) != 0);
}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::begin(const SynchronizationTag & tag,
                              const CommunicationSendRecv & sr) {
  CommunicationPerProcs & comms = this->getCommunications(tag, sr);
  return iterator(this->schemes[sr].begin(), comms.begin(), *this, tag);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::end(const SynchronizationTag & tag,
                            const CommunicationSendRecv & sr) {
  CommunicationPerProcs & comms = this->getCommunications(tag, sr);
  return iterator(this->schemes[sr].end(), comms.end(), *this, tag);
}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAny(const SynchronizationTag & tag,
                                const CommunicationSendRecv & sr) {
  CommunicationPerProcs & comms = this->getCommunications(tag, sr);
  CommunicationPerProcs::iterator it = comms.begin();
  CommunicationPerProcs::iterator end = comms.end();

  std::vector<CommunicationRequest> requests;

  for (; it != end; ++it) {
    requests.push_back(it->second.request);
  }

  UInt req_id = communicator.waitAny(requests);
  if (req_id != UInt(-1)) {
    CommunicationRequest & request = requests[req_id];
    UInt proc =
        sr == _recv ? request.getSource() : request.getDestination();

    return iterator(this->schemes[sr].find(proc),
                    comms.find(proc), *this, tag);
  } else {
    return this->end(tag, sr);
  }
}

/* ---------------------------------------------------------------------- */
template <class Entity>
void
Communications<Entity>::waitAll(const SynchronizationTag & tag,
                                const CommunicationSendRecv & sr) {
  CommunicationPerProcs & comms = this->getCommunications(tag, sr);
  CommunicationPerProcs::iterator it = comms.begin();
  CommunicationPerProcs::iterator end = comms.end();

  std::vector<CommunicationRequest> requests;

  for (; it != end; ++it) {
    requests.push_back(it->second.request);
  }

  communicator.waitAll(requests);
}

template <class Entity>
void Communications<Entity>::incrementPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) {
  ++(pending_communications[sr][tag]);
}

template <class Entity>
void Communications<Entity>::decrementPending(
    const SynchronizationTag & tag, const CommunicationSendRecv & sr) {
  ++(pending_communications[sr][tag]);
}

template <class Entity>
void Communications<Entity>::freeRequests(const SynchronizationTag & tag,
                                          const CommunicationSendRecv & sr) {
  iterator it = this->begin(tag, sr);
  iterator end = this->end(tag, sr);

  for (; it != end; ++it) {
    (*it).freeRequest();
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createScheme(UInt proc,
                                     const CommunicationSendRecv & sr) {
  scheme_iterator it = schemes[sr].find(proc);
  if (it != schemes[sr].end()) {
    AKANTU_CUSTOM_EXCEPTION_INFO(debug::CommunicationException(),
                                 "Communication scheme("
                                     << sr
                                     << ") already created for proc: " << proc);
  }
  return schemes[sr][proc];
}

template <class Entity>
void Communications<Entity>::resetSchemes(const CommunicationSendRecv & sr) {
  scheme_iterator it = this->schemes[sr].begin();
  scheme_iterator end = this->schemes[sr].end();
  for (; it != end; ++it) {
    it->second.resize(0);
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::initializeCommunication(
    const SynchronizationTag & tag, UInt proc, UInt size,
    const CommunicationSendRecv & sr) {
  // accessor that creates if it does not exists
  CommunicationsPerTags & comms = communications[sr];
  CommunicationPerProcs & comms_per_tag = comms[tag];
  Communication & comm = comms_per_tag[proc];
  comm.size = size;
  comm.type = sr;

  if (comm_counter.find(tag) == comm_counter.end()) {
    comm_counter[tag] = 0;
  }

  if (pending_communications[sr].find(tag) ==
      pending_communications[sr].end()) {
    pending_communications[sr][tag] = 0;
  }
}

/* -------------------------------------------------------------------------- */
template <class Entity>
Communications<Entity>::Communications(StaticCommunicator & communicator)
    : communicator(communicator) {}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::begin_send(const SynchronizationTag & tag) {
  return this->begin(tag, _send);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::end_send(const SynchronizationTag & tag) {
  return this->end(tag, _send);
}

/* ---------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::begin_recv(const SynchronizationTag & tag) {
  return this->begin(tag, _recv);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::end_recv(const SynchronizationTag & tag) {
  return this->end(tag, _recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::tag_iterator
Communications<Entity>::begin_tag() {
  return tag_iterator(comm_counter.begin());
}

template <class Entity>
typename Communications<Entity>::tag_iterator
Communications<Entity>::end_tag() {
  return tag_iterator(comm_counter.end());
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::begin_send_scheme() {
  return this->schemes[_send].begin();
}

template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::end_send_scheme() {
  return this->schemes[_send].end();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::begin_send_scheme() const {
  return this->schemes[_send].begin();
}

template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::end_send_scheme() const {
  return this->schemes[_send].end();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::begin_recv_scheme() {
  return this->schemes[_recv].begin();
}

template <class Entity>
typename Communications<Entity>::scheme_iterator
Communications<Entity>::end_recv_scheme() {
  return this->schemes[_recv].end();
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::begin_recv_scheme() const {
  return this->schemes[_recv].begin();
}

template <class Entity>
typename Communications<Entity>::const_scheme_iterator
Communications<Entity>::end_recv_scheme() const {
  return this->schemes[_recv].end();
}

/* ---------------------------------------------------------------------- */
template <class Entity>
bool Communications<Entity>::hasCommunication(
    const SynchronizationTag & tag) const {
  return (communications[_send].find(tag) != communications[_send].end());
}

template <class Entity>
UInt Communications<Entity>::incrementCounter(const SynchronizationTag & tag) {
  std::map<SynchronizationTag, UInt>::iterator it = comm_counter.find(tag);
  if (it == comm_counter.end()) {
    AKANTU_CUSTOM_EXCEPTION_INFO(
        debug::CommunicationException(),
        "No known communications for the tags: " << tag);
  }

  ++(it->second);
  return it->second;
}

template <class Entity>
bool Communications<Entity>::hasPendingRecv(
    const SynchronizationTag & tag) const {
  return this->hasPending(tag, _recv);
}

template <class Entity>
bool Communications<Entity>::hasPendingSend(
    const SynchronizationTag & tag) const {
  return this->hasPending(tag, _send);
}

template <class Entity>
StaticCommunicator & Communications<Entity>::getCommunicator() const {
  return communicator;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAnyRecv(const SynchronizationTag & tag) {
  return this->waitAny(tag, _recv);
}

template <class Entity>
typename Communications<Entity>::iterator
Communications<Entity>::waitAnySend(const SynchronizationTag & tag) {
  return this->waitAny(tag, _send);
}

template <class Entity>
void Communications<Entity>::waitAllRecv(const SynchronizationTag & tag) {
  this->waitAll(tag, _recv);
}

template <class Entity>
void Communications<Entity>::waitAllSend(const SynchronizationTag & tag) {
  this->waitAll(tag, _send);
}

template <class Entity>
void Communications<Entity>::freeSendRequests(const SynchronizationTag & tag) {
  this->freeRequests(tag, _send);
}

template <class Entity>
void Communications<Entity>::freeRecvRequests(const SynchronizationTag & tag) {
  this->freeRequests(tag, _recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createSendScheme(UInt proc) {
  return createScheme(proc, _send);
}

template <class Entity>
typename Communications<Entity>::Scheme &
Communications<Entity>::createRecvScheme(UInt proc) {
  return createScheme(proc, _recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity> void Communications<Entity>::resetSchemes() {
  resetSchemes(_send);
  resetSchemes(_recv);
}

/* -------------------------------------------------------------------------- */
template <class Entity>
const typename Communications<Entity>::Scheme &
Communications<Entity>::getSendScheme(UInt proc) const {
  return this->schemes[_send].find(proc)->seconf;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
const typename Communications<Entity>::Scheme &
Communications<Entity>::getRecvScheme(UInt proc) const {
  return this->schemes[_recv].find(proc)->seconf;
}

/* -------------------------------------------------------------------------- */
template <class Entity>
void Communications<Entity>::initializeSendCommunication(
    const SynchronizationTag & tag, UInt proc, UInt size) {
  initializeCommunication(tag, proc, size, _send);
}

template <class Entity>
void Communications<Entity>::initializeRecvCommunication(
    const SynchronizationTag & tag, UInt proc, UInt size) {
  initializeCommunication(tag, proc, size, _send);
}

} // akantu

#endif /* __AKANTU_COMMUNICATIONS_TMPL_HH__ */
