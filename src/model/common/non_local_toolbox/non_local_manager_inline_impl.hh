/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH_
#define AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void NonLocalManager::registerNeighborhood(const ID & neighborhood,
                                                  const ID & weight_func_id) {

  /// check if neighborhood has already been created
  auto it = neighborhoods.find(neighborhood);
  if (it == neighborhoods.end()) {
    this->createNeighborhood(weight_func_id, neighborhood);
  }
}

/* -------------------------------------------------------------------------- */
inline NonLocalNeighborhoodBase &
NonLocalManager::getNeighborhood(const ID & name) const {
  AKANTU_DEBUG_IN();

  auto it = neighborhoods.find(name);

  AKANTU_DEBUG_ASSERT(it != neighborhoods.end(),
                      "The neighborhood " << name << " is not registered");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

} // namespace akantu

#endif /* AKANTU_NON_LOCAL_MANAGER_INLINE_IMPL_HH_ */
