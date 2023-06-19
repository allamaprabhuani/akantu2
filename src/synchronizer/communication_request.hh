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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_REAL_STATIC_COMMUNICATOR_HH_
#define AKANTU_REAL_STATIC_COMMUNICATOR_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
class DLL_PUBLIC InternalCommunicationRequest {
public:
  InternalCommunicationRequest(Idx source, Idx dest);
  virtual ~InternalCommunicationRequest();

  virtual void printself(std::ostream & stream, int indent = 0) const;

  AKANTU_GET_MACRO(Source, source, Idx);
  AKANTU_GET_MACRO(Destination, destination, Idx);

private:
  Idx source;
  Idx destination;
  Idx id;
  static Int counter;
};

/* -------------------------------------------------------------------------- */
class CommunicationRequest {
public:
  CommunicationRequest(
      std::shared_ptr<InternalCommunicationRequest> request = nullptr)
      : request(std::move(request)) {}

  virtual ~CommunicationRequest() = default;

  virtual void free() { request.reset(); }

  void printself(std::ostream & stream, int indent = 0) const {
    request->printself(stream, indent);
  };

  Idx getSource() const { return request->getSource(); }
  Idx getDestination() const { return request->getDestination(); }

  bool isFreed() const { return request == nullptr; }

  InternalCommunicationRequest & getInternal() { return *request; }

private:
  std::shared_ptr<InternalCommunicationRequest> request;
};

/* -------------------------------------------------------------------------- */
class CommunicationStatus {
public:
  AKANTU_GET_MACRO(Source, source, Int);
  Int size() const { return size_; }
  AKANTU_GET_MACRO(Tag, tag, Int);

  AKANTU_SET_MACRO(Source, source, Int);
  AKANTU_SET_MACRO(Size, size_, Int);
  AKANTU_SET_MACRO(Tag, tag, Int);

private:
  Idx source{0};
  Int size_{0};
  Int tag{0};
};

/* -------------------------------------------------------------------------- */
/// Datatype to pack pairs for MPI_{MIN,MAX}LOC
template <typename T1, typename T2> struct SCMinMaxLoc {
  T1 min_max;
  T2 loc;
};

} // namespace akantu

#endif /* AKANTU_REAL_STATIC_COMMUNICATOR_HH_ */
