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
#include "model.hh"
//#include "neighborhood_max_criterion.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH_
#define AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Int
NeighborhoodMaxCriterion::getNbDataForElements(const Array<Element> & elements,
                                               SynchronizationTag tag) const {
  auto nb_quadrature_points = this->model.getNbIntegrationPoints(elements);
  Int size = 0;

  if (tag == SynchronizationTag::_nh_criterion) {
    size += sizeof(Real) * nb_quadrature_points;
  }

  return size;
}

/* -------------------------------------------------------------------------- */
inline void
NeighborhoodMaxCriterion::packElementData(CommunicationBuffer & buffer,
                                          const Array<Element> & elements,
                                          SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_nh_criterion) {
    NeighborhoodMaxCriterion::packElementalDataHelper(
        criterion, buffer, elements, true, this->model.getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */
inline void
NeighborhoodMaxCriterion::unpackElementData(CommunicationBuffer & buffer,
                                            const Array<Element> & elements,
                                            SynchronizationTag tag) {
  if (tag == SynchronizationTag::_nh_criterion) {
    NeighborhoodMaxCriterion::unpackElementalDataHelper(
        criterion, buffer, elements, true, this->model.getFEEngine());
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_NEIGHBORHOOD_MAX_CRITERION_INLINE_IMPL_HH_ */
