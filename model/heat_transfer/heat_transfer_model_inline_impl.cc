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

