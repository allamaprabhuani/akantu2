/**
 * @file   heat_transfer_model_inline_impl.cc
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @date   Fri Mar  4 17:04:25 2011
 *
 * @brief  Implementation of the inline functions of the HeatTransferModel class
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

/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToPack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  switch(tag) {
  case _gst_htm_temperature: 
  case _gst_htm_capacity: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  
  AKANTU_DEBUG_OUT();
  return size;  
};
/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToUnpack(SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = getFEM().getMesh().getNbNodes();

  switch(tag) {
  case _gst_htm_capacity: 
  case _gst_htm_temperature: {
    size += nb_nodes * sizeof(Real);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  
  AKANTU_DEBUG_OUT();
  return size;  
};
/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
					const UInt index,
					SynchronizationTag tag) const{
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_htm_capacity: 
    buffer << (*capacity_lumped)(index);
    break;
  case _gst_htm_temperature: {
    buffer << (*temperature)(index);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  
  AKANTU_DEBUG_OUT();
};
/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
					  const UInt index,
					  SynchronizationTag tag) const{

  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_htm_capacity: {
    buffer >> (*capacity_lumped)(index);
    break;
  }
  case _gst_htm_temperature: {
    buffer >> (*temperature)(index);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
  
  AKANTU_DEBUG_OUT();
};
/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToPack(const Element & element,
					       SynchronizationTag tag) const {

  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}
/* -------------------------------------------------------------------------- */
inline UInt HeatTransferModel::getNbDataToUnpack(const Element & element,
						 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return size;
}
/* -------------------------------------------------------------------------- */
inline void HeatTransferModel::packData(CommunicationBuffer & buffer,
					const Element & element,
					SynchronizationTag tag) const {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getConnectivity(element.type).values;
  
  switch (tag){ 
  case _gst_htm_capacity: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer << (*capacity_lumped)(offset_conn);
    }
    break;
  }
  case _gst_htm_temperature: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer << (*temperature)(offset_conn);
    }
    break;
  }
  }
}
  /* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
		       const Element & element,
		       SynchronizationTag tag) const {

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getConnectivity(element.type,_ghost).values;

  switch (tag){ 
  case _gst_htm_capacity: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer >> (*capacity_lumped)(offset_conn);
    }
    break;
  }
  case _gst_htm_temperature: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer >> (*temperature)(offset_conn);
    }  
    break;
  }
  }
}
/* -------------------------------------------------------------------------- */


