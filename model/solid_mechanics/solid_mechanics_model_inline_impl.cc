/**
 * @file   solid_mechanics_model_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 29 12:07:04 2010
 *
 * @brief  Implementation of the inline functions of the SolidMechanicsModel class
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
inline Material & SolidMechanicsModel::getMaterial(UInt mat_index) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(mat_index < materials.size(),
		      "The model " << id << " has no material no "<< mat_index);
  AKANTU_DEBUG_OUT();
  return *materials[mat_index];
}



/* -------------------------------------------------------------------------- */
template <typename M>
UInt SolidMechanicsModel::readCustomMaterial(const std::string & filename,
					     const std::string & keyword) {

  MaterialParser parser;
  parser.open(filename);
  std::string key = keyword;
  to_lower(key);
  std::string mat_name = parser.getNextMaterialType();
  while (mat_name != ""){
    if (mat_name == key) break;
    mat_name = parser.getNextMaterialType();
  }
  if (mat_name != key) AKANTU_DEBUG_ERROR("material "
					  << key
					  << " not found in file " << filename);

  std::stringstream sstr_mat; sstr_mat << id << ":" << materials.size() << ":" << key;
  MaterialID mat_id = sstr_mat.str();
  Material * mat = parser.readMaterialObject<M>(*this,mat_id);
  materials.push_back(mat);
  return materials.size();;
}

/* -------------------------------------------------------------------------- */
inline FEM & SolidMechanicsModel::getFEMBoundary(std::string name) {
  return dynamic_cast<FEM &>(getFEMClassBoundary<MyFEMType>(name));
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataToPack(const Element & element,
						 GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifdef AKANTU_DEBUG
  size += spatial_dimension; /// position of the barycenter of the element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element; // mass vector
    break;
  }
  case _gst_smm_for_strain: {
    size += nb_nodes_per_element * spatial_dimension; // displacement

    UInt mat = element_material[element.type]->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    size += nb_nodes_per_element * spatial_dimension * 3; // force, displacement, boundary
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
inline UInt SolidMechanicsModel::getNbDataToUnpack(const Element & element,
						   GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifdef AKANTU_DEBUG
  size += spatial_dimension; /// position of the barycenter of the element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element; // mass vector
    break;
  }
  case _gst_smm_for_strain: {
    size += nb_nodes_per_element * spatial_dimension; // displacement

    UInt mat = ghost_element_material[element.type]->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    size += nb_nodes_per_element * spatial_dimension * 3; // force, displacement, boundary
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
inline void SolidMechanicsModel::packData(Real ** buffer,
					  const Element & element,
					  GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getConnectivity(element.type).values;

#ifdef AKANTU_DEBUG
  getFEM().getMesh().getBarycenter(element.element, element.type, *buffer);
  (*buffer) += spatial_dimension;
#endif

  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      (*buffer)[n] = mass->values[offset_conn];
    }
    *buffer += nb_nodes_per_element;
    break;
  }
  case _gst_smm_for_strain: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;
      memcpy(*buffer, current_position->values + offset_conn, spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;
    }

    UInt mat = element_material[element.type]->values[element.element];
    materials[mat]->packData(buffer, element, tag);
    break;
  }
  case _gst_smm_boundary: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;

      memcpy(*buffer, force->values + offset_conn, spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;

      memcpy(*buffer, velocity->values + offset_conn, spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;

      for (UInt i = 0; i < spatial_dimension; ++i) {
	(*buffer)[i] = boundary->values[offset_conn + i] ? 1.0 : -1.0;
      }
      *buffer += spatial_dimension;
    }
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(Real ** buffer,
					    const Element & element,
					    GhostSynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = getFEM().getMesh().getGhostConnectivity(element.type).values;

#ifdef AKANTU_DEBUG
  Real barycenter[spatial_dimension];
  getFEM().getMesh().getBarycenter(element.element, element.type, barycenter, _ghost);

  Real tolerance = 1e-15;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    if(!(fabs(barycenter[i] - (*buffer)[i]) <= tolerance))
      AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element : "
			 << element
			 << "(barycenter[" << i << "] = " << barycenter[i]
			 << " and (*buffer)[" << i << "] = " << (*buffer)[i] << ")");
  }
  *buffer += spatial_dimension;
#endif

  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      mass->values[offset_conn] = (*buffer)[n];
    }
    *buffer += nb_nodes_per_element;
    break;
  }
  case _gst_smm_for_strain: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;
      memcpy(current_position->values + offset_conn, *buffer,  spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;
    }

    UInt mat = ghost_element_material[element.type]->values[element.element];
    materials[mat]->unpackData(buffer, element, tag);
    break;
  }
  case _gst_smm_boundary: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n] * spatial_dimension;

      memcpy(force->values + offset_conn, *buffer, spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;

      memcpy(velocity->values + offset_conn, *buffer, spatial_dimension * sizeof(Real));
      *buffer += spatial_dimension;

      for (UInt i = 0; i < spatial_dimension; ++i) {
	boundary->values[offset_conn + i] = ((*buffer)[i] > 0) ? true : false;
      }
      *buffer += spatial_dimension;
    }
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}
