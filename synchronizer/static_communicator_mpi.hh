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
#include "real_static_communicator.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class CommunicationRequestMPI : public CommunicationRequest {
public:
  __aka_inline__ CommunicationRequestMPI(UInt source, UInt dest);
  __aka_inline__ ~CommunicationRequestMPI();
  __aka_inline__ MPI_Request * getMPIRequest() { return request; };
private:
  MPI_Request * request;
};


/* -------------------------------------------------------------------------- */
class StaticCommunicatorMPI : public RealStaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  __aka_inline__ StaticCommunicatorMPI(int * argc, char *** argv);

  __aka_inline__ virtual ~StaticCommunicatorMPI();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  template<typename T> __aka_inline__ void send(T * buffer, Int size, Int receiver, Int tag);
  template<typename T> __aka_inline__ void receive(T * buffer, Int size, Int sender, Int tag);

  template<typename T> __aka_inline__ CommunicationRequest * asyncSend(T * buffer, Int size, Int receiver, Int tag);
  template<typename T> __aka_inline__ CommunicationRequest * asyncReceive(T * buffer, Int size, Int sender, Int tag);

  template<typename T> __aka_inline__ void allGather(T * values, Int nb_values);
  template<typename T> __aka_inline__ void allGatherV(T * values, Int * nb_values);

  template<typename T> __aka_inline__ void gather(T * values, Int nb_values, Int root);
  template<typename T> __aka_inline__ void gatherV(T * values, Int * nb_values, Int root);
  template<typename T> __aka_inline__ void broadcast(T * values, Int nb_values, Int root);

  __aka_inline__ bool testRequest(CommunicationRequest * request);

  __aka_inline__ void wait(CommunicationRequest * request);

  __aka_inline__ void waitAll(std::vector<CommunicationRequest *> & requests);

  __aka_inline__ void barrier();

  template<typename T> __aka_inline__ void allReduce(T * values, Int nb_values, const SynchronizerOperation & op);

private:
  template<typename T>
  __aka_inline__ MPI_Datatype getMPIDatatype();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  __aka_inline__ void setMPICommunicator(MPI_Comm comm);
  __aka_inline__ MPI_Comm getMPICommunicator() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  MPI_Comm communicator;

  static MPI_Op synchronizer_operation_to_mpi_op[_so_null + 1];
};


/* -------------------------------------------------------------------------- */
/* __aka_inline__ functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "static_communicator_mpi_inline_impl.cc"
#endif


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_MPI_HH__ */
