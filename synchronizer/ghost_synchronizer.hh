/**
 * @file   ghost_synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 20 17:40:08 2010
 *
 * @brief  Class of ghost synchronisation (PBC or parallel communication)
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GHOST_SYNCHRONIZER_HH__
#define __AKANTU_GHOST_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class Synchronizer;
}

__BEGIN_AKANTU__

class GhostSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  GhostSynchronizer();

  virtual ~GhostSynchronizer();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// synchronize ghosts
  void synchronize(GhostSynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(GhostSynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(GhostSynchronizationTag tag);

  /// register a new communication
  void registerTag(GhostSynchronizationTag tag, const std::string & name);

  /// register a new synchronization
  void registerSynchronizer(Synchronizer & synchronizer);

public:
  virtual UInt getNbDataToPack(const Element & element,
			       GhostSynchronizationTag tag) const = 0;

  virtual UInt getNbDataToUnpack(const Element & element,
				 GhostSynchronizationTag tag) const = 0;

  virtual void packData(Real ** buffer,
			const Element & element,
			GhostSynchronizationTag tag) const = 0;

  virtual void unpackData(Real ** buffer,
			  const Element & element,
			  GhostSynchronizationTag tag) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// list of synchronizers
  std::list<Synchronizer *> synchronizers;

  /// number of tags registered
  UInt nb_ghost_synchronization_tags;

  /// list of registered synchronization
  std::map<GhostSynchronizationTag, std::string> registered_synchronization;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "ghost_synchronizer_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const GhostSynchronizer & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_GHOST_SYNCHRONIZER_HH__ */
