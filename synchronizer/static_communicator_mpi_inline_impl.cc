/**
 * @file   static_communicator_mpi_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep  2 15:10:51 2010
 *
 * @brief  implementation of the mpi communicator
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

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
inline StaticCommunicatorMPI::StaticCommunicatorMPI(int * argc, char *** argv) {
  MPI_Init(argc, argv);
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
inline MPI_Comm StaticCommunicatorMPI::getMPICommunicator() {
  return communicator;
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(UInt * buffer, Int size,
					Int receiver, Int tag) {
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Send(buffer, size, MPI_UNSIGNED, receiver, tag, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::send(Real * buffer, Int size,
					Int receiver, Int tag) {
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Send(buffer, size, MPI_DOUBLE, receiver, tag, communicator);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Send.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(UInt * buffer, Int size,
					   Int sender, Int tag) {
  MPI_Status status;

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Recv(buffer, size, MPI_UNSIGNED, sender, tag, communicator, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::receive(Real * buffer, Int size,
					   Int sender, Int tag) {
  MPI_Status status;

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Recv(buffer, size, MPI_DOUBLE, sender, tag, communicator, &status);

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Recv.");
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest * StaticCommunicatorMPI::asyncSend(UInt * buffer, Int size,
							     Int receiver, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(prank, receiver);

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Isend(buffer, size, MPI_UNSIGNED, receiver, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest * StaticCommunicatorMPI::asyncSend(Real * buffer, Int size,
							     Int receiver, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(prank, receiver);

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Isend(buffer, size, MPI_DOUBLE, receiver, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Isend.");
  return request;
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest * StaticCommunicatorMPI::asyncReceive(UInt * buffer, Int size,
								  Int sender, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(sender, prank);

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Irecv(buffer, size, MPI_UNSIGNED, sender, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Irecv.");
  return request;
}

/* -------------------------------------------------------------------------- */
inline CommunicationRequest * StaticCommunicatorMPI::asyncReceive(Real * buffer, Int size,
								  Int sender, Int tag) {
  CommunicationRequestMPI * request = new CommunicationRequestMPI(sender, prank);

#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Irecv(buffer, size, MPI_DOUBLE, sender, tag, communicator, request->getMPIRequest());

  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Irecv.");
  return request;
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
inline void StaticCommunicatorMPI::allReduce(UInt * values, UInt nb_values,
					     const SynchronizerOperation & op) {
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allreduce(MPI_IN_PLACE, values, nb_values, MPI_UNSIGNED,
		  synchronizer_operation_to_mpi_op[op],
		  communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allreduce.");
}

/* -------------------------------------------------------------------------- */
inline void StaticCommunicatorMPI::allReduce(Real * values, UInt nb_values,
					     const SynchronizerOperation & op) {
#if !defined(AKANTU_NDEBUG)
  int ret =
#endif
    MPI_Allreduce(MPI_IN_PLACE, values, nb_values, MPI_DOUBLE,
		  synchronizer_operation_to_mpi_op[op],
		  communicator);
  AKANTU_DEBUG_ASSERT(ret == MPI_SUCCESS, "Error in MPI_Allreduce.");
}
