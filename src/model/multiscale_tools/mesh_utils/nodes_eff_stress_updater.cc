/**
 * @file   nodes_eff_stress_updater.cc
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @date Tue Feb 8  2022
 * @brief  synchronizer for the effective stresses at the shared nodes
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
/* ------------------------------------------------------------------ */
#include "nodes_eff_stress_updater.hh"
#include "element_synchronizer.hh"
/* ------------------------------------------------------------------- */
#include <numeric>
/* ------------------------------------------------------------------- */

namespace akantu {

void NodesEffStressUpdater::updateMaxEffStressAtNodes() {
  if (mesh.getCommunicator().getNbProc() == 1)
    return;
  this->synchronizer.synchronizeOnce(*this,
                                     SynchronizationTag::_rve_border_nodes);
  this->synchronizer.slaveReductionOnce(*this,
                                        SynchronizationTag::_rve_border_nodes);
}

} // namespace akantu
