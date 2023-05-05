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

} // namespace akantu

#include "igfem_helper.hh"
#include "solid_mechanics_model_igfem.hh"
#include <iostream>

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt MaterialIGFEM::getNbDataForElements(const Array<Element> & elements,
                                                SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension *
           spatial_dimension * sizeof(Real) *
           this->getModel().getNbIntegrationPoints(elements, "IGFEMFEEngine");
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void MaterialIGFEM::packElementData(CommunicationBuffer & buffer,
                                           const Array<Element> & elements,
                                           SynchronizationTag tag) const {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      packElementDataHelper(piola_kirchhoff_2, buffer, elements,
                            "IGFEMFEEngine");
      packElementDataHelper(gradu, buffer, elements, "IGFEMFEEngine");
    }
    packElementDataHelper(stress, buffer, elements, "IGFEMFEEngine");
  }
}

/* -------------------------------------------------------------------------- */
inline void MaterialIGFEM::unpackElementData(CommunicationBuffer & buffer,
                                             const Array<Element> & elements,
                                             SynchronizationTag tag) {
  if (tag == SynchronizationTag::_smm_stress) {
    if (this->isFiniteDeformation()) {
      unpackElementDataHelper(piola_kirchhoff_2, buffer, elements,
                              "IGFEMFEEngine");
      unpackElementDataHelper(gradu, buffer, elements, "IGFEMFEEngine");
    }
    unpackElementDataHelper(stress, buffer, elements, "IGFEMFEEngine");
  }
}
