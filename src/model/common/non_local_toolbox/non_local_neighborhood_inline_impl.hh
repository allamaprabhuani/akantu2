/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include "non_local_neighborhood.hh"
#include "non_local_manager.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_HH_
#define AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_HH_

namespace akantu {
/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline Int NonLocalNeighborhood<WeightFunction>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  Int size = 0;

  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto && variable_id : non_local_variables) {
      size += this->non_local_manager.getNbData(elements, variable_id);
    }
  }

  size += this->weight_function->getNbData(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto && variable_id : non_local_variables) {
      this->non_local_manager.packData(buffer, elements, variable_id);
    }
  }

  this->weight_function->packData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto & variable_id : non_local_variables) {
      this->non_local_manager.unpackData(buffer, elements, variable_id);
    }
  }

  this->weight_function->unpackData(buffer, elements, tag);
}

} // namespace akantu

#endif /* AKANTU_NON_LOCAL_NEIGHBORHOOD_INLINE_IMPL_HH_ */
