/**
 * @file   synchronizer.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 23 13:48:37 2010
 *
 * @brief  interface for communicator and pbc synchronizers
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

#ifndef __AKANTU_SYNCHRONIZER_HH__
#define __AKANTU_SYNCHRONIZER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "data_accessor.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
  class GhostSynchronizer;
}

__BEGIN_AKANTU__

class Synchronizer : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Synchronizer(SynchronizerID id = "synchronizer", MemoryID memory_id = 0);

  virtual ~Synchronizer() { };

  virtual void printself(std::ostream & stream, int indent = 0) const {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// synchronize ghosts
  void synchronize(DataAccessor & data_accessor,SynchronizationTag tag);

  /// asynchronous synchronization of ghosts
  virtual void asynchronousSynchronize(DataAccessor & data_accessor,SynchronizationTag tag) = 0;

  /// wait end of asynchronous synchronization of ghosts
  virtual void waitEndSynchronize(DataAccessor & data_accessor,SynchronizationTag tag) = 0;

  /// compute buffer size for a given tag and data accessor
  virtual void computeBufferSize(DataAccessor & data_accessor, SynchronizationTag tag)=0;

protected:

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


/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const Synchronizer & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
