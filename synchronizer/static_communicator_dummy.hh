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

  virtual void send(__attribute__ ((unused)) UInt * buffer,
		    __attribute__ ((unused)) Int size,
		    __attribute__ ((unused)) Int receiver,
		    __attribute__ ((unused)) Int tag) {};
  virtual void send(__attribute__ ((unused)) Real * buffer,
		    __attribute__ ((unused)) Int size,
		    __attribute__ ((unused)) Int receiver,
		    __attribute__ ((unused)) Int tag) {};

  virtual void receive(__attribute__ ((unused)) UInt * buffer,
		       __attribute__ ((unused)) Int size,
		       __attribute__ ((unused)) Int sender,
		       __attribute__ ((unused)) Int tag) {};
  virtual void receive(__attribute__ ((unused)) Real * buffer,
		       __attribute__ ((unused)) Int size,
		       __attribute__ ((unused)) Int sender,
		       __attribute__ ((unused)) Int tag) {};

  virtual CommunicationRequest * asyncSend(__attribute__ ((unused)) UInt * buffer,
					   __attribute__ ((unused)) Int size,
					   __attribute__ ((unused)) Int receiver,
					   __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncSend(__attribute__ ((unused)) Real * buffer,
					   __attribute__ ((unused)) Int size,
					   __attribute__ ((unused)) Int receiver,
					   __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(__attribute__ ((unused)) UInt * buffer,
					      __attribute__ ((unused)) Int size,
					      __attribute__ ((unused)) Int sender,
					      __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual CommunicationRequest * asyncReceive(__attribute__ ((unused)) Real * buffer,
					      __attribute__ ((unused)) Int size,
					      __attribute__ ((unused)) Int sender,
					      __attribute__ ((unused)) Int tag) {
    return new CommunicationRequest(0, 0);
  };

  virtual bool testRequest(__attribute__ ((unused)) CommunicationRequest * request) { return true; };


  virtual void wait(__attribute__ ((unused)) CommunicationRequest * request) {};

  virtual void waitAll(__attribute__ ((unused)) std::vector<CommunicationRequest *> & requests) {};

  virtual void barrier() {};

  virtual void allReduce(__attribute__ ((unused)) Real * values,
			 __attribute__ ((unused)) UInt nb_values,
			 __attribute__ ((unused)) const SynchronizerOperation & op) {};
  virtual void allReduce(__attribute__ ((unused)) UInt * values,
			 __attribute__ ((unused)) UInt nb_values,
			 __attribute__ ((unused)) const SynchronizerOperation & op) {};

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
