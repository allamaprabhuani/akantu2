/**
 * @file   static_communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:34:09 2010
 *
 * @brief  Class handling the parallel communications
 *
 * @section LICENSE
 *
 * <insert license here>
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
  virtual ~CommunicationRequest() {};
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
  // {
  //   return new CommunicationRequest(); };
  virtual CommunicationRequest * asyncSend(Real * buffer, Int size,
					   Int receiver, Int tag) = 0;
  // {
  //   return new CommunicationRequest();
  // };

  virtual void wait(CommunicationRequest * request) = 0;

  virtual void waitAll(std::vector<CommunicationRequest *> & requests) = 0;

  virtual void freeCommunicationRequest(CommunicationRequest * request);
  virtual void freeCommunicationRequest(std::vector<CommunicationRequest *> & requests);

  virtual void barrier() = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() = 0;
  virtual Int whoAmI() = 0;

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


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
