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
inline const Material & SolidMechanicsModel::getMaterial(UInt mat_index) const {
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

  Parser parser;
  parser.open(filename);
  std::string key = keyword;
  to_lower(key);
  std::string mat_name = parser.getNextSection("material");
  while (mat_name != ""){
    if (mat_name == key) break;
    mat_name = parser.getNextSection("material");
  }
  if (mat_name != key) AKANTU_DEBUG_ERROR("material "
					  << key
					  << " not found in file " << filename);

  std::stringstream sstr_mat; sstr_mat << id << ":" << materials.size() << ":" << key;
  MaterialID mat_id = sstr_mat.str();
  Material * mat = parser.readSection<M>(*this,mat_id);
  materials.push_back(mat);
  return materials.size();;
}

/* -------------------------------------------------------------------------- */
inline FEM & SolidMechanicsModel::getFEMBoundary(std::string name) {
  return dynamic_cast<FEM &>(getFEMClassBoundary<MyFEMType>(name));
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataToPack(const Element & element,
						 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifndef AKANTU_NDEBUG
  size += spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) * spatial_dimension; // mass vector
    break;
  }
  case _gst_smm_for_strain: {
    size += nb_nodes_per_element * spatial_dimension * sizeof(Real); // displacement

    UInt mat = element_material(element.type, _not_ghost)->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
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
						   SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);

#ifndef AKANTU_NDEBUG
  size += spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
#endif

  switch(tag) {
  case _gst_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) * spatial_dimension; // mass vector
    break;
  }
  case _gst_smm_for_strain: {
    size += nb_nodes_per_element * spatial_dimension * sizeof(Real); // displacement

    UInt mat = element_material(element.type, _ghost)->values[element.element];
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
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
inline void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
					  const Element & element,
					  SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _not_ghost;

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = mesh.getConnectivity(element.type, ghost_type).values;

#ifndef AKANTU_NDEBUG
  types::RVector barycenter(spatial_dimension);
  mesh.getBarycenter(element.element, element.type, barycenter.storage(), ghost_type);
  buffer << barycenter;
#endif

  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      Vector<Real>::iterator<types::RVector> it_mass = mass->begin(spatial_dimension);
      buffer << it_mass[offset_conn];
    }
    break;
  }
  case _gst_smm_for_strain: {
    Vector<Real>::iterator<types::RVector> it_disp = displacement->begin(spatial_dimension);
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer << it_disp[offset_conn];
    }

    UInt mat = element_material(element.type, ghost_type)->values[element.element];
    materials[mat]->packData(buffer, element, tag);
    break;
  }
  case _gst_smm_boundary: {
    Vector<Real>::iterator<types::RVector> it_force = force->begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> it_velocity = velocity->begin(spatial_dimension);
    Vector<bool>::iterator<types::Vector<bool> > it_boundary = boundary->begin(spatial_dimension);

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];

      buffer << it_force   [offset_conn];
      buffer << it_velocity[offset_conn];
      buffer << it_boundary[offset_conn];
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
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
					    const Element & element,
					    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  GhostType ghost_type = _ghost;

  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(element.type);
  UInt el_offset  = element.element * nb_nodes_per_element;
  UInt * conn  = mesh.getConnectivity(element.type, ghost_type).values;

#ifndef AKANTU_NDEBUG
  types::RVector barycenter_loc(spatial_dimension);
  mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(), ghost_type);

  types::RVector barycenter(spatial_dimension);
  buffer >> barycenter;
  Real tolerance = 1e-15;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    if(!(std::abs(barycenter(i) - barycenter_loc(i)) <= tolerance))
      AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element : "
			 << element
			 << "(barycenter[" << i << "] = " << barycenter_loc(i)
			 << " and buffer[" << i << "] = " << barycenter(i) << ")");
  }
#endif

  switch(tag) {
  case _gst_smm_mass: {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      Vector<Real>::iterator<types::RVector> it_mass = mass->begin(spatial_dimension);
      buffer >> it_mass[offset_conn];
    }
    break;
  }
  case _gst_smm_for_strain: {
    Vector<Real>::iterator<types::RVector> it_disp = displacement->begin(spatial_dimension);
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];
      buffer >> it_disp[offset_conn];
    }

    UInt mat = element_material(element.type, ghost_type)->values[element.element];
    materials[mat]->unpackData(buffer, element, tag);
    break;
  }
  case _gst_smm_boundary: {
    Vector<Real>::iterator<types::RVector> it_force = force->begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> it_velocity = velocity->begin(spatial_dimension);
    Vector<bool>::iterator<types::Vector<bool> > it_boundary = boundary->begin(spatial_dimension);

    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt offset_conn = conn[el_offset + n];

      buffer >> it_force   [offset_conn];
      buffer >> it_velocity[offset_conn];
      buffer >> it_boundary[offset_conn];
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
inline UInt SolidMechanicsModel::getNbDataToPack(SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += nb_nodes * sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_mass: {
    size += nb_nodes * sizeof(Real) * spatial_dimension;
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
inline UInt SolidMechanicsModel::getNbDataToUnpack(SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += nb_nodes * sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_mass: {
    size += nb_nodes * sizeof(Real) * spatial_dimension;
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
inline void SolidMechanicsModel::packData(CommunicationBuffer & buffer,
					  const UInt index,
					  SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_smm_uv: {
    for (UInt d = 0; d < spatial_dimension; ++d) {
      buffer << (*displacement)(index,d);
      buffer << (*velocity)(index,d);
    }
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("pack mass of node " << index << " which is " << (*mass)(index,0));
    buffer << (*mass)(index,0);
    buffer << (*mass)(index,1);
    buffer << (*mass)(index,2);
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
					    const UInt index,
					    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_smm_uv: {
    for (UInt d = 0; d < spatial_dimension; ++d) {
      buffer >> (*displacement)(index,d);
      buffer >> (*velocity)(index,d);
    }
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("mass of node " << index << " was " << (*mass)(index,0));
    buffer >> (*mass)(index,0);
    buffer >> (*mass)(index,1);
    buffer >> (*mass)(index,2);
    AKANTU_DEBUG_INFO("mass of node " << index << " is now " << (*mass)(index,0));
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

