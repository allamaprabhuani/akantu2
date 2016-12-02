/**
 * @file   synchronizer_registry.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 09 2014
 *
 * @brief  Registry of synchronizers
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
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <map>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SYNCHRONIZER_REGISTRY_HH__
#define __AKANTU_SYNCHRONIZER_REGISTRY_HH__

namespace akantu {
class DataAccessorBase;
class Synchronizer;
}

namespace akantu {

class SynchronizerRegistry {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SynchronizerRegistry();
  virtual ~SynchronizerRegistry();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// synchronize operation
  void synchronize(SynchronizationTag tag);

  /// asynchronous synchronization
  void asynchronousSynchronize(SynchronizationTag tag);

  /// wait end of asynchronous synchronization
  void waitEndSynchronize(SynchronizationTag tag);

  /// register a new synchronization
  void registerSynchronizer(Synchronizer & synchronizer,
                            SynchronizationTag tag);

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;
protected:
  /// Register a different data accessor.
  void registerDataAccessor(DataAccessorBase & data_accessor);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::multimap<SynchronizationTag, Synchronizer *> Tag2Sync;
  /// list of registered synchronization
  Tag2Sync synchronizers;

  /// data accessor that will permit to do the pack/unpack things
  DataAccessorBase * data_accessor;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

// #include "synchronizer_registry_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const SynchronizerRegistry & _this) {
  _this.printself(stream);
  return stream;
}

} // akantu

#endif /* __AKANTU_SYNCHRONIZER_REGISTRY_HH__ */
