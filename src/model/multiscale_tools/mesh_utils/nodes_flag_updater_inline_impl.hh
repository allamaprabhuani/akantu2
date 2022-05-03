/**
 * @file   nodes_flag_updater.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  inline implementation of the synchronizer for the crack-employed
 * nodes flag
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
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "nodes_flag_updater.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt NodesFlagUpdater::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_rve_border_nodes) {
    size +=
        Mesh::getNbNodesPerElementList(elements) * sizeof(bool) + sizeof(int);
  }
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NodesFlagUpdater::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_rve_border_nodes) {
    return;
  }

  int prank = mesh.getCommunicator().whoAmI();
  buffer << prank;
}

/* -------------------------------------------------------------------------- */
inline void NodesFlagUpdater::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_rve_border_nodes) {
    return;
  }

  MeshAccessor mesh_accessor(mesh);
  auto dim = mesh.getSpatialDimension();
  auto & nodes_pos = mesh_accessor.getNodes();
  auto pos_it = nodes_pos.begin(dim);

  int proc;
  buffer >> proc;

  // looping through all the slave elements and preventing insertion at their
  // nodes; the problem could be that a ghost element could have a non-ghost
  // node
  for (const auto & element : elements) {
    // get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    // loop on all connectivity nodes
    for (auto node : current_conn) {
      prevent_insertion(node) = true;
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_NODES_FLAG_UPDATER_INLINE_IMPL_CC__ */
