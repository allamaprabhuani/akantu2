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

  virtual void send(UInt * buffer, Int size, Int receiver, Int tag) { };
  virtual void send(Real * buffer, Int size, Int receiver, Int tag) { };

  virtual void receive(UInt * buffer, Int size, Int sender, Int tag) { };
  virtual void receive(Real * buffer, Int size, Int sender, Int tag) { };

  virtual CommunicationRequest asyncSend(UInt * buffer, Int size,
					 Int receiver, Int tag) {
    return CommunicationRequest(); };
  virtual CommunicationRequest asyncSend(Real * buffer, Int size,
					 Int receiver, Int tag) {
    return CommunicationRequest(); };

  virtual void wait(CommunicationRequest & request) {};

  virtual void waitAll(std::vector<CommunicationRequest> & requests) {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  static StaticCommunicator * getStaticCommunicator(CommunicatorType type = _communicator_mpi);

  static StaticCommunicator * getStaticCommunicator(int * argc, char *** argv, CommunicatorType type = _communicator_mpi);

  static bool isInstantiated() { return is_instantiated; };

  virtual Int getNbProc() { return 1; };

  virtual Int whoAmI() { return 0; };
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

//#include "static_communicator_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
