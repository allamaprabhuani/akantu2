/**
 * @file   nodes_eff_stress_updater_inline_impl.hh
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  inline implementation for the synchronizer for the effective stresses
 * at the shared nodes
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
/* ------------------------------------------------------------------- */
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "nodes_eff_stress_updater.hh"
/* ------------------------------------------------------------------- */

#ifndef __AKANTU_NODES_EFF_STRESS_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_NODES_EFF_STRESS_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* ------------------------------------------------------------------- */
inline UInt
NodesEffStressUpdater::getNbData(const Array<Element> & elements,
                                 const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_border_nodes) {
    size +=
        Mesh::getNbNodesPerElementList(elements) * sizeof(Real) + sizeof(int);
  }
  return size;
}

/* ------------------------------------------------------------------- */
inline void
NodesEffStressUpdater::packData(CommunicationBuffer & buffer,
                                const Array<Element> & elements,
                                const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_border_nodes) {
    return;
  }

  int prank = mesh.getCommunicator().whoAmI();
  buffer << prank;

  for (const auto & element : elements) {
    /// get element connectivity
    const Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      buffer << nodes_eff_stress(node);
    }
  }
}

/* ------------------------------------------------------------------ */
inline void NodesEffStressUpdater::unpackData(CommunicationBuffer & buffer,
                                              const Array<Element> & elements,
                                              const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_border_nodes) {
    return;
  }

  MeshAccessor mesh_accessor(mesh);

  int proc;
  buffer >> proc;

  for (const auto & element : elements) {
    /// get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      Real received_eff_stress;
      buffer >> received_eff_stress;
      auto master_or_slave =
          (mesh.isMasterNode(node) or mesh.isSlaveNode(node));
      if (master_or_slave and (nodes_eff_stress(node) < received_eff_stress)) {
        nodes_eff_stress(node) = 0;
      }
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_NODES_EFF_STRESS_UPDATER_INLINE_IMPL_CC__ */
