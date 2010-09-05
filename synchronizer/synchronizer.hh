/**
 * @file   synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 23 13:48:37 2010
 *
 * @brief  interface for communicator and pbc synchronizers
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_HH__
#define __AKANTU_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class Synchronizer : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Synchronizer(SynchronizerID id = "synchronizer", MemoryID memory_id = 0);

  virtual ~Synchronizer() { };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// function to print the contain of the class
  //  virtual void printself(std::ostream & stream, int indent = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// id of the synchronizer
  SynchronizerID id;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "synchronizer_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const Synchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
