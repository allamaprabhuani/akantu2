/**
 * @file   static_communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:34:09 2010
 *
 * @brief  Class handling the parallel communications
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_HH__
#define __AKANTU_STATIC_COMMUNICATOR_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class CommunicationRequest {
public:
  CommunicationRequest(UInt source, UInt dest);
  virtual ~CommunicationRequest();

  virtual void printself(std::ostream & stream, int indent = 0) const;

  AKANTU_GET_MACRO(Source, source, UInt);
  AKANTU_GET_MACRO(Destination, destination, UInt);
private:
  UInt source;
  UInt destination;
  UInt id;
  static UInt counter;
};


class StaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  StaticCommunicator() { };

public:
  virtual ~StaticCommunicator() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void send(UInt * buffer, Int size, Int receiver, Int tag) = 0;
  virtual void send(Real * buffer, Int size, Int receiver, Int tag) = 0;

  virtual void receive(UInt * buffer, Int size, Int sender, Int tag) = 0;
  virtual void receive(Real * buffer, Int size, Int sender, Int tag) = 0;

  virtual CommunicationRequest * asyncSend(UInt * buffer, Int size,
					   Int receiver, Int tag) = 0;
  virtual CommunicationRequest * asyncSend(Real * buffer, Int size,
					   Int receiver, Int tag) = 0;

  virtual CommunicationRequest * asyncReceive(UInt * buffer, Int size,
					      Int sender, Int tag) = 0;
  virtual CommunicationRequest * asyncReceive(Real * buffer, Int size,
					      Int sender, Int tag) = 0;

  virtual bool testRequest(CommunicationRequest * request) = 0;

  virtual void wait(CommunicationRequest * request) = 0;

  virtual void waitAll(std::vector<CommunicationRequest *> & requests) = 0;

  virtual void freeCommunicationRequest(CommunicationRequest * request);
  virtual void freeCommunicationRequest(std::vector<CommunicationRequest *> & requests);

  virtual void barrier() = 0;

  virtual void allReduce(Real * values, Int nb_values, const SynchronizerOperation & op) = 0;
  virtual void allReduce(UInt * values, Int nb_values, const SynchronizerOperation & op) = 0;

  virtual void gather(Real * values, Int nb_values, Int root) = 0;
  virtual void gather(UInt * values, Int nb_values, Int root) = 0;
  virtual void gather(Int * values, Int nb_values, Int root) = 0;

  virtual void gatherv(Real * values, Int * nb_values, Int root) = 0;
  virtual void gatherv(UInt * values, Int * nb_values, Int root) = 0;
  virtual void gatherv(Int * values, Int * nb_values, Int root) = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() const = 0;
  virtual Int whoAmI() const = 0;

  static StaticCommunicator * getStaticCommunicator(CommunicatorType type = _communicator_mpi);

  static StaticCommunicator * getStaticCommunicator(int * argc, char *** argv, CommunicatorType type = _communicator_mpi);

  static bool isInstantiated() { return is_instantiated; };
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  static bool is_instantiated;

  static StaticCommunicator * static_communicator;

protected:
  Int prank;

  Int psize;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "static_communicator_inline_impl.cc"

/* -------------------------------------------------------------------------- */
/* Inline Functions VectorBase                                                */
/* -------------------------------------------------------------------------- */
inline std::ostream & operator<<(std::ostream & stream, const CommunicationRequest & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
