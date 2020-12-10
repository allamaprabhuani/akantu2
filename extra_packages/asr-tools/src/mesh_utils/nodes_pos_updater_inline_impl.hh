/**
 * @file   global_ids_updater_inline_impl.hh
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
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "nodes_pos_updater.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NODES_POS_UPDATER_INLINE_IMPL_CC__
#define __AKANTU_NODES_POS_UPDATER_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt NodesPosUpdater::getNbData(const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  UInt size = 0;
  if (tag == SynchronizationTag::_asr) {
    auto dim = mesh.getSpatialDimension();
    size += Mesh::getNbNodesPerElementList(elements) * dim * sizeof(Real) +
            Mesh::getNbNodesPerElementList(elements) * sizeof(bool) +
            sizeof(int);
  }
  return size;
}

/* -------------------------------------------------------------------------- */
inline void NodesPosUpdater::packData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      const SynchronizationTag & tag) const {
  if (tag != SynchronizationTag::_asr)
    return;

  int prank = mesh.getCommunicator().whoAmI();
  buffer << prank;

  auto dim = mesh.getSpatialDimension();
  auto & positions = mesh.getNodes();
  auto pos_it = positions.begin(dim);

  for (auto & element : elements) {
    /// get element connectivity
    const Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      Vector<Real> coord(pos_it[node]);
      buffer << coord;
      buffer << modified_pos(node);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void NodesPosUpdater::unpackData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        const SynchronizationTag & tag) {
  if (tag != SynchronizationTag::_asr)
    return;

  MeshAccessor mesh_accessor(mesh);
  auto dim = mesh.getSpatialDimension();
  auto & nodes_pos = mesh_accessor.getNodes();
  auto pos_it = nodes_pos.begin(dim);

  int proc;
  buffer >> proc;

  for (auto & element : elements) {
    /// get element connectivity
    Vector<UInt> current_conn =
        const_cast<const Mesh &>(mesh).getConnectivity(element);

    /// loop on all connectivity nodes
    for (auto node : current_conn) {
      Vector<Real> current_pos(pos_it[node]);
      auto current_modif_flag = modified_pos(node);
      Vector<Real> unpacked_pos(dim);
      bool unpacked_mod_flag;
      buffer >> unpacked_pos;
      buffer >> unpacked_mod_flag;
      if (this->reduce) {
        /// if slave is modified -> take position of slave; if both slave and
        /// master are modified take the average; if master is modified or
        /// noneof them are modified -> do nothing
        if (unpacked_mod_flag and not current_modif_flag) {
          current_pos = unpacked_pos;
        } else if (unpacked_mod_flag and current_modif_flag) {
          current_pos += unpacked_pos;
          current_pos /= 2;
        }
      } else {
        current_pos = unpacked_pos;
      }
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_NODES_POS_UPDATER_INLINE_IMPL_CC__ */
