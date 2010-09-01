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

__BEGIN_AKANTU__

class StaticCommunicator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  StaticCommunicator() { };

public:
  virtual ~StaticCommunicator() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  StaticCommunicator * getStaticCommunicator();

  static bool isInstantiated() const { return isInstantiated; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  static isInstantiated;

  static StaticCommunicator * static_communicator;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "static_communicator_inline_impl.cc"


__END_AKANTU__

#endif /* __AKANTU_STATIC_COMMUNICATOR_HH__ */
