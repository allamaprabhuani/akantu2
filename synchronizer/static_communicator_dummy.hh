/**
 * @file   static_communicator_dummy.hh
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

#ifndef __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__
#define __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class StaticCommunicatorDummy : public StaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  virtual ~StaticCommunicatorDummy() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void send(UInt * buffer, Int size, Int receiver, Int tag) {};
  virtual void send(Real * buffer, Int size, Int receiver, Int tag) {};

  virtual void receive(UInt * buffer, Int size, Int sender, Int tag) {};
  virtual void receive(Real * buffer, Int size, Int sender, Int tag) {};

  virtual CommunicationRequest * asyncSend(UInt * buffer, Int size,
					   Int receiver, Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncSend(Real * buffer, Int size,
					   Int receiver, Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(UInt * buffer, Int size,
					      Int sender, Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(Real * buffer, Int size,
					      Int sender, Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual bool testRequest(CommunicationRequest * request) { return true; };


  virtual void wait(CommunicationRequest * request) {};

  virtual void waitAll(std::vector<CommunicationRequest *> & requests) {};

  virtual void barrier() {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual Int getNbProc() { return 1; };
  virtual Int whoAmI() { return 0; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
};

__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_DUMMY_HH__ */
