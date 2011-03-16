/**
 * @file   static_communicator_mpi.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep  2 19:59:58 2010
 *
 * @brief  class handling parallel communication trough MPI
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STATIC_COMMUNICATOR_MPI_HH__
#define __AKANTU_STATIC_COMMUNICATOR_MPI_HH__

/* -------------------------------------------------------------------------- */
#include <mpi.h>
/* -------------------------------------------------------------------------- */
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class CommunicationRequestMPI : public CommunicationRequest {
public:
  inline CommunicationRequestMPI(UInt source, UInt dest);
  inline ~CommunicationRequestMPI();
  inline MPI_Request * getMPIRequest() { return request; };
private:
  MPI_Request * request;
};


/* -------------------------------------------------------------------------- */
class StaticCommunicatorMPI : public virtual StaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  inline StaticCommunicatorMPI(int * argc, char *** argv);

  inline virtual ~StaticCommunicatorMPI();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  inline void send(UInt * buffer, Int size, Int receiver, Int tag) {
    _send(buffer, size, receiver, tag);
  };
  inline void send(Int * buffer, Int size, Int receiver, Int tag) {
    _send(buffer, size, receiver, tag);
  };
  inline void send(Real * buffer, Int size, Int receiver, Int tag) {
    _send(buffer, size, receiver, tag);
  };

  inline void receive(UInt * buffer, Int size, Int sender, Int tag) {
    _receive(buffer, size, sender, tag);
  };
  inline void receive(Int * buffer, Int size, Int sender, Int tag) {
    _receive(buffer, size, sender, tag);
  };
  inline void receive(Real * buffer, Int size, Int sender, Int tag) {
    _receive(buffer, size, sender, tag);
  };

  inline CommunicationRequest * asyncSend(UInt * buffer, Int size, Int receiver, Int tag) {
    return _asyncSend(buffer, size, receiver, tag);
  };
  inline CommunicationRequest * asyncSend(Int * buffer, Int size, Int receiver, Int tag) {
    return _asyncSend(buffer, size, receiver, tag);
  };
  inline CommunicationRequest * asyncSend(Real * buffer, Int size, Int receiver, Int tag) {
    return _asyncSend(buffer, size, receiver, tag);
  };

  inline CommunicationRequest * asyncReceive(UInt * buffer, Int size,
					      Int sender, Int tag) {
    return _asyncReceive(buffer, size, sender, tag);
  };
  inline CommunicationRequest * asyncReceive(Real * buffer, Int size,
					      Int sender, Int tag) {
    return _asyncReceive(buffer, size, sender, tag);
  };

  inline bool testRequest(CommunicationRequest * request);

  inline void wait(CommunicationRequest * request);

  inline void waitAll(std::vector<CommunicationRequest *> & requests);

  inline void barrier();

  inline void allReduce(Real * values, Int nb_values, const SynchronizerOperation & op);
  inline void allReduce(UInt * values, Int nb_values, const SynchronizerOperation & op);

  inline void gather(Real * values, Int nb_values, Int root = 0) { _gather(values, nb_values, root); };
  inline void gather(UInt * values, Int nb_values, Int root = 0) { _gather(values, nb_values, root); };
  inline void gather(Int  * values, Int nb_values, Int root = 0) { _gather(values, nb_values, root); };

  inline void gatherv(Real * values, Int * nb_values, Int root = 0) { _gatherv(values, nb_values, root); };
  inline void gatherv(UInt * values, Int * nb_values, Int root = 0) { _gatherv(values, nb_values, root); };
  inline void gatherv(Int  * values, Int * nb_values, Int root = 0) { _gatherv(values, nb_values, root); };

  inline void broadcast(Real * values, Int nb_values, Int root = 0) { _broadcast(values, nb_values, root); };
  inline void broadcast(UInt * values, Int nb_values, Int root = 0) { _broadcast(values, nb_values, root); };
  inline void broadcast(Int  * values, Int nb_values, Int root = 0) { _broadcast(values, nb_values, root); };


private:
  template<typename T>
  inline MPI_Datatype getMPIDatatype();

  template<typename T>
  inline void _send(T * buffer, Int size, Int receiver, Int tag);

  template<typename T>
  inline void _receive(T * buffer, Int size, Int sender, Int tag);

  template<typename T>
  inline CommunicationRequest * _asyncSend(T * buffer, Int size, Int receiver, Int tag);

  template<typename T>
  inline CommunicationRequest * _asyncReceive(T * buffer, Int size, Int sender, Int tag);

  template<typename T>
  inline void _gather(T * values, Int nb_values, Int root);

  template<typename T>
  inline void _gatherv(T * values, Int * nb_values, Int root);

  template<typename T>
  inline void _broadcast(T * values, Int nb_values, Int root);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline void setMPICommunicator(MPI_Comm comm);
  inline MPI_Comm getMPICommunicator() const;

  inline Int getNbProc() const { return psize; };
  inline Int whoAmI() const { return prank; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  MPI_Comm communicator;

  static MPI_Op synchronizer_operation_to_mpi_op[_so_null + 1];
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "static_communicator_mpi_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_MPI_HH__ */
