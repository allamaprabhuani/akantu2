/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "global_ids_updater.hh"
#include "element_synchronizer.hh"
#include "mesh_accessor.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

namespace akantu {

std::pair<Int, Int> GlobalIdsUpdater::updateGlobalIDs() {
  // if (mesh.getCommunicator().getNbProc() == 1) {
  //   return {local_nb_new_nodes, };
  // }

  auto && [total_nb_new_nodes, total_nb_new_elements] = this->updateGlobalIDsLocally();

  if (mesh.isDistributed()) {
    this->synchronizeGlobalIDs();
  }
  return {total_nb_new_nodes, total_nb_new_elements};
}

std::pair<Int, Int> GlobalIdsUpdater::updateGlobalIDsLocally() {
  return mesh.updateOffsets();
}

void GlobalIdsUpdater::synchronizeGlobalIDs() {
  this->reduce = true;
  this->synchronizer->slaveReductionOnce(*this,
                                         SynchronizationTag::_giu_global_conn);

#ifndef AKANTU_NDEBUG
  for (auto node : nodes_flags) {
    auto node_flag = mesh.getNodeFlag(node.first);
    if (node_flag != NodeFlag::_pure_ghost) {
      continue;
    }
    auto n = 0U;

    for (auto & pair : node.second) {
      if (std::get<1>(pair) == NodeFlag::_pure_ghost) {
        ++n;
      }
    }

    if (n == node.second.size()) {
      AKANTU_DEBUG_WARNING(
          "The node " << n << "is ghost on all the neighboring processors");
    }
  }
#endif

  this->reduce = false;
  this->synchronizer->synchronizeOnce(*this,
                                      SynchronizationTag::_giu_global_conn);
}

} // namespace akantu
