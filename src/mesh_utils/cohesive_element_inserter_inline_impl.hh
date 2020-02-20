/**
 * @file   cohesive_element_inserter_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Oct 13 2017
 * @date last modification: Wed Nov 11 2020
 *
 * @brief  Cohesive element inserter inline functions
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
//#include "cohesive_element_inserter.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_HH_
#define AKANTU_COHESIVE_ELEMENT_INSERTER_INLINE_IMPL_HH_

namespace akantu {

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
