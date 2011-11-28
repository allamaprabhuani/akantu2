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
					  SynchronizationTag tag) {

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

#ifndef AKANTU_NDEBUG
  size += spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
  size += spatial_dimension * nb_nodes_per_element * sizeof(Real); /// position of the nodes of the element
#endif

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case _gst_htm_gradient_temperature: {
    size += spatial_dimension * sizeof(Real); // temperature gradient
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    size += spatial_dimension * nb_nodes_per_element * sizeof(Real); // shape derivatives
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

#ifndef AKANTU_NDEBUG
  size += spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
  size += spatial_dimension * nb_nodes_per_element * sizeof(Real); /// position of the nodes of the element
#endif

  switch(tag) {
  case _gst_htm_capacity: {
    size += nb_nodes_per_element * sizeof(Real); // capacity vector
    break;
  }
  case _gst_htm_temperature: {
    size += nb_nodes_per_element * sizeof(Real); // temperature
    break;
  }
  case _gst_htm_gradient_temperature: {
    size += spatial_dimension * sizeof(Real); // temperature gradient
    size += nb_nodes_per_element * sizeof(Real); // nodal temperatures
    size += spatial_dimension * nb_nodes_per_element * sizeof(Real); // shape derivatives
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
  GhostType ghost_type = _not_ghost;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getConnectivity(element.type, ghost_type).values;

#ifndef AKANTU_NDEBUG
  types::RVector barycenter(spatial_dimension);
  getFEM().getMesh().getBarycenter(element.element, element.type, barycenter.storage(), ghost_type);
  buffer << barycenter;
  Real * nodes = getFEM().getMesh().getNodes().values;
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt offset_conn = conn[el_offset + n];
    for (UInt s = 0; s < spatial_dimension; ++s) {
    buffer << nodes[spatial_dimension*offset_conn+s];
    }
  }
#endif

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
  case _gst_htm_gradient_temperature: {
    Vector<Real>::const_iterator<types::RVector> it_gtemp =
      temperature_gradient(element.type, ghost_type).begin(spatial_dimension);

    buffer << it_gtemp[element.element];

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer << (*temperature)(offset_conn);
    }

    Vector<Real>::const_iterator<types::Matrix> it_shaped =
      getFEM().getShapesDerivatives(element.type, ghost_type).begin(nb_nodes_per_element,spatial_dimension);
    buffer << it_shaped[element.element];
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}
  /* -------------------------------------------------------------------------- */
inline void HeatTransferModel::unpackData(CommunicationBuffer & buffer,
		       const Element & element,
		       SynchronizationTag tag) {
  GhostType ghost_type = _ghost;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getConnectivity(element.type, ghost_type).values;

#ifndef AKANTU_NDEBUG
  types::RVector barycenter_loc(spatial_dimension);
  getFEM().getMesh().getBarycenter(element.element, element.type, barycenter_loc.storage(), ghost_type);

  types::RVector barycenter(spatial_dimension);
  buffer >> barycenter;
  Real tolerance = 1e-15;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    if(!(std::abs(barycenter(i) - barycenter_loc(i)) <= tolerance))
      AKANTU_EXCEPTION("Unpacking an unknown value for the element : "
			     << element
			     << "(barycenter[" << i << "] = " << barycenter_loc(i)
			     << " and buffer[" << i << "] = " << barycenter(i) << ")");
  }

  types::RVector coords(spatial_dimension);
  Real * nodes = getFEM().getMesh().getNodes().values;
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    buffer >> coords;
    UInt offset_conn = conn[el_offset + n];
    Real tolerance = 1e-15;
    Real * coords_local = nodes+spatial_dimension*offset_conn;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs(coords(i) - coords_local[i]) <= tolerance))
	AKANTU_EXCEPTION("Unpacking to wrong node for the element : "
			 << element
			 << "(coords[" << i << "] = " << coords_local[i]
			 << " and buffer[" << i << "] = " << coords(i) << ")");
    }
  }
#endif

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
      Real tbuffer;
      buffer >> tbuffer;
#ifndef AKANTU_NDEBUG
      // if (!getFEM().getMesh().isPureGhostNode(offset_conn)){
      // 	if (std::abs(tbuffer - (*temperature)(offset_conn)) > 1e-15){
      // 	  AKANTU_EXCEPTION(std::scientific << std::setprecision(20)
      // 			   << "local node is impacted with a different value computed from a distant proc"
      // 			   << " => divergence of trajectory detected " << std::endl
      // 			   << tbuffer << " != " << (*temperature)(offset_conn)
      // 			   << " diff is " << tbuffer - (*temperature)(offset_conn));
      // 	}
      // }
#endif
    (*temperature)(offset_conn) = tbuffer;
    }
    break;
  }
  case _gst_htm_gradient_temperature: {
    Vector<Real>::iterator<types::RVector> it_gtemp =
      temperature_gradient(element.type, ghost_type).begin(spatial_dimension);
    types::RVector gtemp(spatial_dimension);
    types::RVector temp_nodes(nb_nodes_per_element);
    types::Matrix shaped(nb_nodes_per_element,spatial_dimension);

    buffer >> gtemp;
    buffer >> temp_nodes;
    buffer >> shaped;

    //    Real tolerance = 1e-15;
    if (!Math::are_vector_equal(spatial_dimension,gtemp.storage(),it_gtemp[element.element].storage())){
      Real dist = Math::distance_3d(gtemp.storage(), it_gtemp[element.element].storage());
      debug::debugger.getOutputStream().precision(20);
      std::stringstream temperatures_str;
      temperatures_str.precision(20);
      temperatures_str << std::scientific << "temperatures are ";
      for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	UInt offset_conn = conn[el_offset + n];
	temperatures_str << (*temperature)(offset_conn) << " ";
      }
      Vector<Real>::iterator<types::Matrix> it_shaped =
	const_cast<Vector<Real> &>(getFEM().getShapesDerivatives(element.type, ghost_type))
	.begin(nb_nodes_per_element,spatial_dimension);


      AKANTU_EXCEPTION("packed gradient do not match for element " << element.element << std::endl
		       << "buffer is " << gtemp << " local is " << it_gtemp[element.element]
		       << " dist is " << dist << std::endl
		       << temperatures_str.str() << std::endl
		       << std::scientific << std::setprecision(20)
		       << " distant temperatures " << temp_nodes
		       << "temperature gradient size " << temperature_gradient(element.type, ghost_type).getSize()
		       << " number of ghost elements " << getFEM().getMesh().getNbElement(element.type,_ghost)
		       << std::scientific << std::setprecision(20)
		       << " shaped " << shaped
		       << std::scientific << std::setprecision(20)
		       << " local shaped " << it_shaped[element.element]);
    }
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }
}
/* -------------------------------------------------------------------------- */
