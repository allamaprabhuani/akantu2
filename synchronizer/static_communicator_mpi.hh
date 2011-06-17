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
  inline CommunicationRequestMPI(UInt source, UInt dest);
  inline ~CommunicationRequestMPI();
  inline MPI_Request * getMPIRequest() { return request; };
private:
  MPI_Request * request;
};


/* -------------------------------------------------------------------------- */
class StaticCommunicatorMPI : public RealStaticCommunicator {
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

  template<typename T> inline void send(T * buffer, Int size, Int receiver, Int tag);
  template<typename T> inline void receive(T * buffer, Int size, Int sender, Int tag);

  template<typename T> inline CommunicationRequest * asyncSend(T * buffer, Int size, Int receiver, Int tag);
  template<typename T> inline CommunicationRequest * asyncReceive(T * buffer, Int size, Int sender, Int tag);

  template<typename T> inline void allGather(T * values, Int nb_values);
  template<typename T> inline void allGatherV(T * values, Int * nb_values);

  template<typename T> inline void gather(T * values, Int nb_values, Int root);
  template<typename T> inline void gatherV(T * values, Int * nb_values, Int root);
  template<typename T> inline void broadcast(T * values, Int nb_values, Int root);

  inline bool testRequest(CommunicationRequest * request);

  inline void wait(CommunicationRequest * request);

  inline void waitAll(std::vector<CommunicationRequest *> & requests);

  inline void barrier();

  template<typename T> inline void allReduce(T * values, Int nb_values, const SynchronizerOperation & op);

private:
  template<typename T>
  inline MPI_Datatype getMPIDatatype();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline void setMPICommunicator(MPI_Comm comm);
  inline MPI_Comm getMPICommunicator() const;

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
