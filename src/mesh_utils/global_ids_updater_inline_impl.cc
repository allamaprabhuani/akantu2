/**
 * @file   global_ids_updater_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 02 2015
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Implementation of the inline functions of GlobalIdsUpdater
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
#include "global_ids_updater.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt GlobalIdsUpdater::getNbData(const Array<Element> & elements,
                                        const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == _gst_giu_global_conn)
    size += Mesh::getNbNodesPerElementList(elements) * sizeof(UInt);

  return size;
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  if (tag != _gst_giu_global_conn)
    return;

  auto & global_nodes_ids = mesh.getGlobalNodesIds();

  for (auto & element : elements) {
    /// get element connectivity
    const Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      UInt index = -1;
      if ((this->reduce and mesh.isLocalOrMasterNode(node)) or
          (not this->reduce and not mesh.isPureGhostNode(node))) {
        index = global_nodes_ids(node);
      }
      buffer << index;
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void GlobalIdsUpdater::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  if (tag != _gst_giu_global_conn)
    return;

  auto & global_nodes_ids = mesh.getGlobalNodesIds();

  for (auto & element : elements) {
    /// get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      UInt index;
      buffer >> index;

      if (index == UInt(-1))
        continue;

      if(mesh.isSlaveNode(node))
        global_nodes_ids(node) = index;

      // if (global_nodes_ids(node) == UInt(-1) and mesh.isSlaveNode(node)) {
      //   global_nodes_ids(node) = index;
      // } else {
      //   if (global_nodes_ids(node) != index) {
      //     //change_event.getList().push_back(node);
      //     //change_event.getFormerGlobalID().push_back(global_nodes_ids(node));
      //     global_nodes_ids(node) = index;
      //   }
      // }
    }
  }
}

} // akantu

#endif /* __AKANTU_GLOBAL_IDS_UPDATER_INLINE_IMPL_CC__ */
