/**
 * @file   ghost_synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Aug 20 17:40:08 2010
 *
 * @brief  Class of ghost synchronisation (PBC or parallel communication)
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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

  /* ------------------------------------------------------------------------ */
  /// synchronize ghosts
  void synchronize(GhostSynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  void asynchronousSynchronize(GhostSynchronizationTag tag);

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronize(GhostSynchronizationTag tag);

  /// reduce a value (essentially for communications synchronizer)
  void allReduce(Real * values, const SynchronizerOperation & op, UInt nb_values = 1);

  /* ------------------------------------------------------------------------ */
  /// register a new communication
  void registerTag(GhostSynchronizationTag tag, const std::string & name);

  /// register a new synchronization
  void registerSynchronizer(Synchronizer & synchronizer);

public:
  /**
   * @brief get  the number of  data to send  for a given akantu::Element  and a
   * given akantu::GhostSynchronizationTag
   */
  virtual UInt getNbDataToPack(const Element & element,
			       GhostSynchronizationTag tag) const = 0;


  /**
   * @brief get the number of data  to receive for a given akantu::Element and a
   * given akantu::GhostSynchronizationTag
   */
  virtual UInt getNbDataToUnpack(const Element & element,
				 GhostSynchronizationTag tag) const = 0;

  /**
   * @brief   pack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::GhostSynchronizationTag
   */
  virtual void packData(Real ** buffer,
			const Element & element,
			GhostSynchronizationTag tag) const = 0;

  /**
   * @brief   unpack  the   data  for   a  given   akantu::Element  and   a  given
   * akantu::GhostSynchronizationTag
   */
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
