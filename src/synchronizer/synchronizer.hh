/**
 * @file   synchronizer.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Thu Dec 10 2015
 *
 * @brief  Common interface for synchronizers
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "static_communicator.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Base class for synchronizers                                               */
/* -------------------------------------------------------------------------- */
class Synchronizer : protected Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Synchronizer(
      const ID & id = "synchronizer", MemoryID memory_id = 0,
      StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator());

  virtual ~Synchronizer(){};

  virtual void printself(__attribute__((unused)) std::ostream & stream,
                         __attribute__((unused)) int indent = 0) const {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// synchronize ghosts
  template<class DataAccessor>
  void synchronize(DataAccessor & data_accessor,
                   const SynchronizationTag & tag);

  /// asynchronous synchronization of ghosts
  template<class DataAccessor>
  void asynchronousSynchronize(DataAccessor & data_accessor,
                               const SynchronizationTag & tag);

  /// wait end of asynchronous synchronization of ghosts
  template<class DataAccessor>
  void waitEndSynchronize(DataAccessor & data_accessor,
                                  const SynchronizationTag & tag);

  /// compute buffer size for a given tag and data accessor
  template<class DataAccessor>
  void computeBufferSize(DataAccessor & data_accessor,
                         const SynchronizationTag & tag);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Communicator, communicator, StaticCommunicator &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id of the synchronizer
  SynchronizerID id;

  /// hashed version of the id
  int hash_id;

  /// message counter per tag
  std::map<SynchronizationTag, UInt> tag_counter;

  /// the static memory instance
  StaticCommunicator & communicator;

  /// nb processors in the communicator
  UInt nb_proc;

  /// rank in the communicator
  UInt rank;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const Synchronizer & _this) {
  _this.printself(stream);
  return stream;
}

} // akantu

#include "synchronizer_tmpl.hh"

#endif /* __AKANTU_SYNCHRONIZER_HH__ */
