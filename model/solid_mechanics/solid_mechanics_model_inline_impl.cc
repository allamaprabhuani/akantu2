/**
 * @file   solid_mechanics_model_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 12:07:04 2010
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel class
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
		      "The model " << id << " has no material no "<< mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataToPack(const Element & element,
						 GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifdef AKANTU_DEBUG
  size += spatial_dimension; /// position of the barycenter of the P1 element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element; // mass vector
    break;
  }
  case _gst_smm_residual: {
    UInt mat = element_material[element.type]->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
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
inline UInt SolidMechanicsModel::getNbDataToUnpack(const Element & element,
						   GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifdef AKANTU_DEBUG
  size += spatial_dimension; /// position of the barycenter of the P1 element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element; // mass vector
    break;
  }
  case _gst_smm_residual: {
    UInt mat = ghost_element_material[element.type]->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
  return 0;
};

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::packData(Real ** buffer,
					  const Element & element,
					  GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;

#ifdef AKANTU_DEBUG
  Real * coord = fem->getMesh().getNodes().values;
  UInt * conn  = fem->getMesh().getConnectivity(element.type).values;
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);

  memset(*buffer, 0, spatial_dimension * sizeof(Real));
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt offset_conn = conn[el_offset + n] * spatial_dimension;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      *buffer[i] += coord[offset_conn + i] / (double) spatial_dimension;
    }
  }
  *buffer += spatial_dimension;
#endif


  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;
      **buffer = mass->values[offset_conn];
      *buffer++;
    }
    break;
  }
  case _gst_smm_residual: {
    UInt mat = ghost_element_material[element.type]->values[element.element];
    materials[mat]->packData(buffer, element, tag);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
};

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(Real ** buffer,
					    const Element & element,
					    GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;

#ifdef AKANTU_DEBUG
  Real * coord = fem->getMesh().getNodes().values;
  UInt * conn  = fem->getMesh().getConnectivity(element.type).values;
  UInt spatial_dimension = Mesh::getSpatialDimension(element.type);

  Real barycenter[spatial_dimension];
  memset(barycenter, 0, spatial_dimension * sizeof(Real));
  for (UInt n = 0; n < nb_nodes_per_element; ++n) {
    UInt offset_conn = conn[el_offset + n] * spatial_dimension;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      barycenter[i] += coord[offset_conn + i] / (double) spatial_dimension;
    }
  }

  Real tolerance = 1e-15;

  for (UInt i = 0; i < spatial_dimension; ++i) {
    if(abs(barycenter[i] - *buffer[i]) <= tolerance)
      AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element : "
			 << element);
  }
  *buffer += spatial_dimension;
#endif

  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;
      mass->values[offset_conn] = **buffer;
      *buffer++;
    }
    break;
  }
  case _gst_smm_residual: {
    UInt mat = ghost_element_material[element.type]->values[element.element];
    materials[mat]->unpackData(buffer, element, tag);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
};
