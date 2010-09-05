/**
 * @file   static_communicator_mpi.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Sep  2 19:59:58 2010
 *
 * @brief  class handling parallel communication trough MPI
 *
 * @section LICENSE
 *
 * <insert license here>
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
  MPI_Request & getMPIRequest() { return request; };
private:
  MPI_Request request;
};


/* -------------------------------------------------------------------------- */
class StaticCommunicatorMPI : public StaticCommunicator {
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

  inline void send(UInt * buffer, Int size, Int receiver, Int tag);
  inline void send(Real * buffer, Int size, Int receiver, Int tag);

  inline void receive(UInt * buffer, Int size, Int sender, Int tag);
  inline void receive(Real * buffer, Int size, Int sender, Int tag);

  inline CommunicationRequest asyncSend(UInt * buffer, Int size, Int receiver, Int tag);
  inline CommunicationRequest asyncSend(Real * buffer, Int size, Int receiver, Int tag);

  inline void wait(CommunicationRequest & request);

  inline void waitAll(std::vector<CommunicationRequest> & requests);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline Int getNbProc() { return psize; };

  inline Int whoAmI() { return prank; };

  inline void setMPIComm(MPI_Comm comm);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  MPI_Comm communicator;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "static_communicator_mpi_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_MPI_HH__ */
