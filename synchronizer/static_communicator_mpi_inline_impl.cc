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
  MPI_Comm_rank(communicator, &prank);
  MPI_Comm_size(communicator, &psize);
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(UInt * buffer, Int size,
					Int receiver, Int tag) {
  Int ret = MPI_Send(buffer, size, MPI_UNSIGNED, receiver, tag, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(Real * buffer, Int size,
					Int receiver, Int tag) {
  Int ret = MPI_Send(buffer, size, MPI_DOUBLE, receiver, tag, communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(UInt * buffer, Int size,
					   Int sender, Int tag) {
  MPI_Status status;
  Int ret = MPI_Recv(buffer, size, MPI_UNSIGNED, sender, tag, communicator, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(Real * buffer, Int size,
					   Int sender, Int tag) {
  MPI_Status status;
  Int ret = MPI_Recv(buffer, size, MPI_DOUBLE, sender, tag, communicator, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest StaticCommunicatorMPI::asyncSend(UInt * buffer, Int size,
							     Int receiver, Int tag) {
  CommunicationRequestMPI request;
  Int ret = MPI_Isend(buffer, size, MPI_UNSIGNED, receiver, tag, communicator, &(request.getMPIRequest()));
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
};

/* -------------------------------------------------------------------------- */
inline CommunicationRequest StaticCommunicatorMPI::asyncSend(Real * buffer, Int size,
							     Int receiver, Int tag) {
  CommunicationRequestMPI request;
  Int ret = MPI_Isend(buffer, size, MPI_DOUBLE, receiver, tag, communicator, &(request.getMPIRequest()));
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
};

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::wait(CommunicationRequest & request) {
  MPI_Status status;
  CommunicationRequestMPI & req_mpi = static_cast<CommunicationRequestMPI &>(request);
  MPI_Request * req = &(req_mpi.getMPIRequest());
  Int ret = MPI_Wait(req, &status);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_WAIT.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::waitAll(std::vector<CommunicationRequest> & requests) {
  MPI_Status status;
  std::vector<CommunicationRequest>::iterator it;
  for(it = requests.begin(); it != requests.end(); ++it) {
    MPI_Request * req = &(static_cast<CommunicationRequestMPI &>(*it).getMPIRequest());
    Int ret = MPI_Wait(req, &status);
    AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_WAIT.");
  }
}
