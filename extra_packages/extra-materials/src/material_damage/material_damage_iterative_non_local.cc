/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_damage_iterative_non_local.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <Int spatial_dimension>
void MaterialDamageIterativeNonLocal<
    spatial_dimension>::computeNonLocalStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if (ghost_type == _not_ghost)
    this->norm_max_equivalent_stress = 0;

  MaterialDamageIterativeNonLocalParent::computeNonLocalStresses(ghost_type);

  /// find global Gauss point with highest stress
  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(this->norm_max_equivalent_stress, SynchronizerOperation::_max);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(damage_iterative_non_local,
                     MaterialDamageIterativeNonLocal);

} // namespace akantu
