/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
//#include "cohesive_element_inserter.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_HH_
#define AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::addPhysicalSurface(const ID & surface_name) {
  physical_surfaces.insert(surface_name);
}

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::addPhysicalVolume(const ID & surface_name) {
  physical_zones.insert(surface_name);
}

/* -------------------------------------------------------------------------- */
inline Int
CohesiveElementInserter::getNbData(const Array<Element> & elements,
                                   const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  Int size = 0;

  if (tag == SynchronizationTag::_ce_groups) {
    size = elements.size() * sizeof(bool);
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::packData(CommunicationBuffer & buffer,
                                  const Array<Element> & elements,
                                  const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();
  if (tag == SynchronizationTag::_ce_groups) {
    for (const auto & el : elements) {
      const bool & data = insertion_facets(el);
      buffer << data;
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void
CohesiveElementInserter::unpackData(CommunicationBuffer & buffer,
                                    const Array<Element> & elements,
                                    const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();

  if (tag == SynchronizationTag::_ce_groups) {
    for (const auto & el : elements) {
      bool & data = insertion_facets(el);
      buffer >> data;
    }
  }
  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_HH_ */
