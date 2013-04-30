/**
 * @file   static_communicator_mpi_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sun Sep 05 16:01:51 2010
 *
 * @brief  implementation of the mpi communicator
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



__END_AKANTU__
#include <iostream>
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
inline CommunicationRequestMPI::CommunicationRequestMPI(UInt source, UInt dest) :
  CommunicationRequest(source, dest) {
  request = new MPI_Request;
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequestMPI::~CommunicationRequestMPI() {
  delete request;
}

/* -------------------------------------------------------------------------- */
inline StaticCommunicatorMPI::StaticCommunicatorMPI(int & argc, char ** & argv) : 
  RealStaticCommunicator(argc, argv) {
  MPI_Init(&argc, &argv);
  setMPICommunicator(MPI_COMM_WORLD);
}

/* -------------------------------------------------------------------------- */
inline StaticCommunicatorMPI::~StaticCommunicatorMPI() {
  MPI_Finalize();
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::setMPICommunicator(MPI_Comm comm) {
  communicator = comm;
  MPI_Comm_rank(communicator, &prank);
  MPI_Comm_size(communicator, &psize);
}

/* -------------------------------------------------------------------------- */
inline MPI_Comm StaticCommunicatorMPI::getMPICommunicator() const {
  return communicator;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::send(T * buffer, Int size,
					Int receiver, Int tag) {
  MPI_Datatype type = getMPIDatatype<T>();
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Send(buffer, size, type, receiver, tag, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::receive(T * buffer, Int size,
					   Int sender, Int tag) {
  MPI_Status status;
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Recv(buffer, size, type, sender, tag, communicator, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationRequest * StaticCommunicatorMPI::asyncSend(T * buffer, Int size,
							       Int receiver, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(prank, receiver);
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Isend(buffer, size, type, receiver, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline CommunicationRequest * StaticCommunicatorMPI::asyncReceive(T * buffer, Int size,
								  Int sender, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(sender, prank);
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Irecv(buffer, size, type, sender, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Irecv.");
  return request;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::probe(Int sender, Int tag,
                                         CommunicationStatus & status) {
  MPI_Status mpi_status;
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Probe(sender, tag, communicator, &mpi_status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Probe.");

  MPI_Datatype type = getMPIDatatype<T>();
  int count;
  MPI_Get_count(&mpi_status, type, &count);

  status.setSource(mpi_status.MPI_SOURCE);
  status.setTag(mpi_status.MPI_TAG);
  status.setSize(count);
}

/* -------------------------------------------------------------------------- */
inline bool StaticCommunicatorMPI::testRequest(CommunicationRequest * request) {
  MPI_Status status;
  int flag;
  CommunicationRequestMPI * req_mpi = static_cast<CommunicationRequestMPI *>(request);
  MPI_Request * req = req_mpi->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif

    MPI_Test(req, &flag, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Test.");
  return (flag != 0);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::wait(CommunicationRequest * request) {
  MPI_Status status;
  CommunicationRequestMPI * req_mpi = static_cast<CommunicationRequestMPI *>(request);
  MPI_Request * req = req_mpi->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Wait(req, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::waitAll(std::vector<CommunicationRequest *> & requests) {
  MPI_Status status;
  std::vector<CommunicationRequest *>::iterator it;
  for(it = requests.begin(); it != requests.end(); ++it) {
    MPI_Request * req = static_cast<CommunicationRequestMPI *>(*it)->getMPIRequest();

#if !defined(AKANTU_NDEBUG)
    int ret =
#endif
      MPI_Wait(req, &status);

    AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Wait.");
  }
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::barrier() {
  MPI_Barrier(communicator);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::allReduce(T * values, Int nb_values,
					     const SynchronizerOperation & op) {
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allreduce(MPI_IN_PLACE, values, nb_values, type,
		  synchronizer_operation_to_mpi_op[op],
		  communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allreduce.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::allGather(T * values, Int nb_values) {
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allgather(MPI_IN_PLACE, nb_values, type, values, nb_values, type, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allgather.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::allGatherV(T * values, Int * nb_values) {
  Int * displs = new Int[psize];
  displs[0] = 0;
  for (Int i = 1; i < psize; ++i) {
    displs[i] = displs[i-1] + nb_values[i-1];
  }

  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allgatherv(MPI_IN_PLACE, *nb_values, type, values, nb_values, displs, type, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  delete [] displs;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::gather(T * values, Int nb_values, Int root) {
  T * send_buf = NULL, * recv_buf = NULL;
  if(prank == root) {
    send_buf = (T *) MPI_IN_PLACE;
    recv_buf = values;
  } else {
    send_buf = values;
  }

  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Gather(send_buf, nb_values, type, recv_buf, nb_values, type, root, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::gatherV(T * values, Int * nb_values, Int root) {
  Int * displs = NULL;
  if(prank == root) {
    displs = new Int[psize];
    displs[0] = 0;
    for (Int i = 1; i < psize; ++i) {
      displs[i] = displs[i-1] + nb_values[i-1];
    }
  }

  T * send_buf = NULL, * recv_buf = NULL;
  if(prank == root) {
    send_buf = (T *) MPI_IN_PLACE;
    recv_buf = values;
  } else send_buf = values;

  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Gatherv(send_buf, *nb_values, type, recv_buf, nb_values, displs, type, root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");

  if(prank == root) {
    delete [] displs;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void StaticCommunicatorMPI::broadcast(T * values, Int nb_values, Int root) {
  MPI_Datatype type = getMPIDatatype<T>();

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Bcast(values, nb_values, type, root, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Gather.");
}

/* -------------------------------------------------------------------------- */
// template<typename T>
// inline MPI_Datatype StaticCommunicatorMPI::getMPIDatatype() {
//   return MPI_DATATYPE_NULL;
// }

template<>
inline MPI_Datatype StaticCommunicatorMPI::getMPIDatatype<char>() {
  return MPI_CHAR;
}

template<>
inline MPI_Datatype StaticCommunicatorMPI::getMPIDatatype<Real>() {
  return MPI_DOUBLE;
}

template<>
inline MPI_Datatype StaticCommunicatorMPI::getMPIDatatype<UInt>() {
  return MPI_UNSIGNED;
}

template<>
inline MPI_Datatype StaticCommunicatorMPI::getMPIDatatype<Int>() {
  return MPI_INT;
}

/* -------------------------------------------------------------------------- */

