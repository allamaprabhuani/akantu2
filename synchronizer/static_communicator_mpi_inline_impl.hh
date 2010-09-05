/**
 * @file   static_communicator_mpi_inline_impl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep  2 15:10:51 2010
 *
 * @brief  implementation of the mpi communicator
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline StaticCommunicatorMPI::StaticCommunicatorMPI(int * argc, char *** argv) {
  MPI_Init(argc, argv);
  setMPIComm(MPI_COMM_WORLD);
};

/* -------------------------------------------------------------------------- */
inline StaticCommunicatorMPI::~StaticCommunicatorMPI() {
  MPI_Finalize();
};

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::setMPIComm(MPI_Comm comm) {
  communicator = comm;
  MPI_Comm_rank(comm, &prank);
  MPI_Comm_size(comm, &psize);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(UInt * buffer, UInt size,
					UInt receiver, UInt tag) {
  int ret = MPI_Send(buffer, size, MPI_UNSIGNED, receiver, tag, comm);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(Real * buffer, UInt size,
					UInt receiver, UInt tag) {
  int ret = MPI_Send(buffer, size, MPI_DOUBLE, receiver, tag, comm);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(UInt * buffer, UInt size,
					   UInt sender, UInt tag) {
  MPI_Status status;
  int ret = MPI_Recv(buffer, size, MPI_UNSIGNED, sender, tag, comm, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(Real * buffer, UInt size,
					   UInt sender, UInt tag) {
  MPI_Status status;
  int ret = MPI_Recv(buffer, size, MPI_DOUBLE, sender, tag, comm, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest StaticCommunicatorMPI::asyncSend(UInt * buffer, UInt size, UInt receiver, UInt tag) {
  CommunicationRequestMPI request;
  int ret = MPI_Isend(buffer, size, MPI_UNSIGNED, receiver, tag, comm, &(request.getMPIRequest()));
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
};

/* -------------------------------------------------------------------------- */
inline CommunicationRequest StaticCommunicatorMPI::asyncSend(Real * buffer, UInt size, UInt receiver, UInt tag) {
  CommunicationRequestMPI request;
  int ret = MPI_Isend(buffer, size, MPI_DOUBLE, receiver, tag, comm, &(request.getMPIRequest()));
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
};

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::wait(const CommunicationRequest & request) {
  MPI_Status status;
  int ret = MPI_Wait(&(request.getMPIRequest()), &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_WAIT.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::waitAll(const std::vector<CommunicationRequest> & requests) {
  MPI_Status status;
  std::vector<CommunicationRequest>::const_iterator it;
  for(it = requests.begin(); it != requests.end(); ++it) {
    int ret = MPI_Wait(&(*it->getMPIRequest()), &status);
    AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_WAIT.");
  }
}
