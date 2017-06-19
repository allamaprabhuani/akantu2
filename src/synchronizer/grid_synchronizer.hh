/**
 * @file   grid_synchronizer.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Synchronizer based on spatial grid
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
#include "element_synchronizer.hh"
#include "synchronizer_registry.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GRID_SYNCHRONIZER_HH__
#define __AKANTU_GRID_SYNCHRONIZER_HH__

namespace akantu {

class Mesh;
template <class T> class SpatialGrid;

class GridSynchronizer : public ElementSynchronizer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  GridSynchronizer(Mesh & mesh, const ID & id = "grid_synchronizer",
                   MemoryID memory_id = 0,
                   const bool register_to_event_manager = true);

public:
  virtual ~GridSynchronizer(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
  *Create the Grid Synchronizer:
  *Compute intersection and send info to neighbours that will be stored in
  *ghosts elements
  */
  template <class E>
  static GridSynchronizer *
  createGridSynchronizer(Mesh & mesh, const SpatialGrid<E> & grid,
                         SynchronizerID id = "grid_synchronizer",
                         SynchronizerRegistry * synch_registry = NULL,
                         const std::set<SynchronizationTag> & tags_to_register =
                             std::set<SynchronizationTag>(),
                         MemoryID memory_id = 0,
                         const bool register_to_event_manager = true);

protected:
  /// Define the tags that will be used in the send and receive instructions
  enum CommTags {
    SIZE_TAG = 0,
    DATA_TAG = 1,
    ASK_NODES_TAG = 2,
    SEND_NODES_TAG = 3
  };

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // akantu

#endif /* __AKANTU_GRID_SYNCHRONIZER_HH__ */
