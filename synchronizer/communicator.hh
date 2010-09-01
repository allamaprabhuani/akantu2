/**
 * @file   communicator.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Aug 19 15:28:35 2010
 *
 * @brief  wrapper to the static communicator
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMUNICATOR_HH__
#define __AKANTU_COMMUNICATOR_HH__

/* -------------------------------------------------------------------------- */
#include "static_communicator.hh"
#include "synchronizer.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Communicator : Synchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Communicator();
  virtual ~Communicator();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:


protected:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the static memory instance
  StaticCommunicator * static_communicator;

  /// size of data to send to each processor by communication tag
  std::map< GhostSynchronizationTag, Vector<UInt> > size_to_send;

  /// size of data to receive form each processor by communication tag
  std::map< GhostSynchronizationTag, Vector<UInt> > size_to_receive;

  Vector<Real> * send_buffer;

  Vector<Real> * receive_buffer;

  /// list of real element to send ordered by type then by receiver processors
  Vector<UInt> *** element_to_send;

  /// list of ghost element to receive ordered by type then by sender processors
  Vector<UInt> *** element_to_recieve;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "communicator_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_COMMUNICATOR_HH__ */
