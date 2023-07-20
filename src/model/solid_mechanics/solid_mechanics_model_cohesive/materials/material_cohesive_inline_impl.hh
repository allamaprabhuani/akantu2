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
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
inline Int MaterialCohesive::addFacet(const Element & element) {
  auto & f_filter = facet_filter(element.type, element.ghost_type);
  f_filter.push_back(element.element);
  return f_filter.size() - 1;
}

/* -------------------------------------------------------------------------- */
template <ElementType type>
void MaterialCohesive::computeNormal(const Array<Real> & /*position*/,
                                     Array<Real> & /*normal*/,
                                     GhostType /*ghost_type*/) {}

/* -------------------------------------------------------------------------- */
inline Int MaterialCohesive::getNbData(const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {

  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    return 2 * spatial_dimension * Int(sizeof(Real)) *
           this->getModel().getNbIntegrationPoints(elements,
                                                   "CohesiveFEEngine");
  }
  case SynchronizationTag::_smmc_damage: {
    return Int(sizeof(Real)) * this->getModel().getNbIntegrationPoints(
                                   elements, "CohesiveFEEngine");
  }
  default: {
  }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::packData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) const {
  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    packInternalFieldHelper(tractions, buffer, elements);
    packInternalFieldHelper(contact_tractions, buffer, elements);
    break;
  }
  case SynchronizationTag::_smmc_damage:
    packInternalFieldHelper(damage, buffer, elements);
    break;
  default: {
  }
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialCohesive::unpackData(CommunicationBuffer & buffer,
                                         const Array<Element> & elements,
                                         const SynchronizationTag & tag) {
  switch (tag) {
  case SynchronizationTag::_smm_stress: {
    unpackInternalFieldHelper(tractions, buffer, elements);
    unpackInternalFieldHelper(contact_tractions, buffer, elements);
    break;
  }
  case SynchronizationTag::_smmc_damage:
    unpackInternalFieldHelper(damage, buffer, elements);
    break;
  default: {
  }
  }
}
} // namespace akantu
