/**
 * @file   global_ids_updater.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 * @date last modification: Fri Dec 08 2017
 *
 * @brief  Class that updates the global ids of new nodes that are
 * inserted in the mesh. The functions in this class must be called
 * after updating the node types
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_NODES_POS_UPDATER_HH__
#define __AKANTU_NODES_POS_UPDATER_HH__

/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "mesh.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
class ElementSynchronizer;
} // namespace akantu

namespace akantu {

class NodesPosUpdater : public DataAccessor<Element> {
public:
  NodesPosUpdater(Mesh & mesh, ElementSynchronizer & synchronizer,
                  const Array<bool> modified_pos)
      : mesh(mesh), synchronizer(synchronizer), modified_pos(modified_pos) {
    AKANTU_DEBUG_ASSERT(mesh.getNbNodes() == modified_pos.size(),
                        "Array modified_pos does not have the same size as "
                        "the number of nodes in the mesh");
  }

  void updateNodesPos();

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// Reference to the mesh to update
  Mesh & mesh;

  /// distributed synchronizer to communicate nodes positions
  ElementSynchronizer & synchronizer;

  /// Tells if a reduction is taking place or not
  bool reduce{false};

  /// which tells which nodes were modified
  Array<bool> modified_pos;
};

} // namespace akantu

#include "nodes_pos_updater_inline_impl.hh"

#endif /* __AKANTU_NODES_POS_UPDATER_HH__ */
