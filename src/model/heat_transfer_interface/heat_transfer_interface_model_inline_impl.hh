/**
 * @file   heat_transfer_interface_model_inline_impl.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Wed May 24 2023
 * @date last modification: Wed May 24 2023
 *
 * @brief  Implementation of the inline functions of the
 * HeatTransferInterfaceModel class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_HEAT_TRANSFER_INTERFACE_MODEL_INLINE_IMPL_HH_
#define AKANTU_HEAT_TRANSFER_INTERFACE_MODEL_INLINE_IMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
inline UInt
HeatTransferInterfaceModel::getNbData(const Array<Element> & elements,
                                      const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  if (elements.empty()) {
    return 0;
  }
  if (elements(0).kind() == _ek_regular) {
    size = HeatTransferModel::getNbData(elements, tag);
  } else if (elements(0).kind() == _ek_cohesive) {
    UInt nb_nodes_per_element = 0;
    Array<Element>::const_iterator<Element> it = elements.begin();
    Array<Element>::const_iterator<Element> end = elements.end();
    for (; it != end; ++it) {
      const Element & el = *it;
      nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
    }
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      size += nb_nodes_per_element * sizeof(Real); // temperature
      break;
    }
    case SynchronizationTag::_htm_temperature_on_qpoints: {
      size +=
          getNbIntegrationPoints(elements, "InterfacesFEEngine") * sizeof(Real);
      break;
    }
    case SynchronizationTag::_htm_gradient_temperature: {
      // temperature gradient
      size += getNbIntegrationPoints(elements, "InterfacesFEEngine") *
              spatial_dimension * sizeof(Real);
      // normal openings
      size +=
          getNbIntegrationPoints(elements, "InterfacesFEEngine") * sizeof(Real);
      // nodal temperatures
      size += nb_nodes_per_element * sizeof(Real);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
  AKANTU_DEBUG_OUT();
  return size;
}

/* --------------------------------------------------------------------------
 */
inline void
HeatTransferInterfaceModel::packData(CommunicationBuffer & buffer,
                                     const Array<Element> & elements,
                                     const SynchronizationTag & tag) const {
  if (elements.empty()) {
    return;
  }
  if (elements(0).kind() == _ek_regular) {
    HeatTransferModel::packData(buffer, elements, tag);
  } else if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      packNodalDataHelper(*temperature, buffer, elements, mesh);
      break;
    }
    case SynchronizationTag::_htm_temperature_on_qpoints: {
      packElementalDataHelper(temperature_on_qpoints_coh, buffer, elements,
                              true, getFEEngine("InterfacesFEEngine"));
      break;
    }
    case SynchronizationTag::_htm_gradient_temperature: {
      packElementalDataHelper(temperature_gradient, buffer, elements, true,
                              getFEEngine("InterfacesFEEngine"));
      packElementalDataHelper(opening_on_qpoints, buffer, elements, true,
                              getFEEngine("InterfacesFEEngine"));
      packNodalDataHelper(*temperature, buffer, elements, mesh);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
}

/* --------------------------------------------------------------------------
 */
inline void
HeatTransferInterfaceModel::unpackData(CommunicationBuffer & buffer,
                                       const Array<Element> & elements,
                                       const SynchronizationTag & tag) {
  if (elements.empty()) {
    return;
  }
  if (elements(0).kind() == _ek_regular) {
    HeatTransferModel::unpackData(buffer, elements, tag);
  } else if (elements(0).kind() == _ek_cohesive) {
    switch (tag) {
    case SynchronizationTag::_htm_temperature: {
      unpackNodalDataHelper(*temperature, buffer, elements, mesh);
      break;
    }
    case SynchronizationTag::_htm_temperature_on_qpoints: {
      unpackElementalDataHelper(temperature_on_qpoints_coh, buffer, elements,
                                true, getFEEngine("InterfacesFEEngine"));
      break;
    }
    case SynchronizationTag::_htm_gradient_temperature: {
      unpackElementalDataHelper(temperature_gradient, buffer, elements, true,
                                getFEEngine("InterfacesFEEngine"));
      unpackElementalDataHelper(opening_on_qpoints, buffer, elements, true,
                                getFEEngine("InterfacesFEEngine"));
      unpackNodalDataHelper(*temperature, buffer, elements, mesh);
      break;
    }
    default: {
      AKANTU_ERROR("Unknown ghost synchronization tag : " << tag);
    }
    }
  }
}

} // namespace akantu

#endif /* AKANTU_HEAT_TRANSFER_MODEL_INLINE_IMPL_HH_ */
