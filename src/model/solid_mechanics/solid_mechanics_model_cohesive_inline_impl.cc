/**
 * @file   solid_mechanics_model_cohesive_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Jan 18 09:40:01 2013
 *
 * @brief  Solid mechanics model inline implementation
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

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelCohesive::splitElementByKind(const Vector<Element> & elements,
							    Vector<Element> & elements_regular,
							    Vector<Element> & elements_cohesive) const {
  for (UInt el = 0; el < elements.getSize(); ++el) {
    if (elements(el).kind == _ek_regular)
      elements_regular.push_back(elements(el));
    else if (elements(el).kind == _ek_cohesive)
      elements_cohesive.push_back(elements(el));
    else
      AKANTU_DEBUG_ERROR("Unrecognized element kind in spliElementByKind.");
  }  
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModelCohesive::getNbDataForElements(const Vector<Element> & elements,
							      SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  Vector<Element> elements_regular;
  Vector<Element> elements_cohesive;

  splitElementByKind(elements, elements_regular, elements_cohesive);

  UInt size_regular = SolidMechanicsModel::getNbDataForElements(elements_regular, tag);

  UInt size = 0;

#ifndef AKANTU_NDEBUG
  size += elements_cohesive.getSize() * spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
#endif

  UInt nb_nodes_per_element = 0;

  Vector<Element>::iterator<Element> it  = elements_cohesive.begin();
  Vector<Element>::iterator<Element> end = elements_cohesive.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch(tag) {
  case _gst_material_id: {
    size += elements_cohesive.getSize() * 2 * sizeof(UInt);
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  default: {  }
  }

  if(tag != _gst_material_id) {
    Vector<Element> * elements_per_mat = new Vector<Element>[materials.size()];
    this->splitElementByMaterial(elements_cohesive, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      size += materials[i]->getNbDataForElements(elements_per_mat[i], tag);
    }
    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
  return size + size_regular;
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelCohesive::packElementData(CommunicationBuffer & buffer,
							 const Vector<Element> & elements,
							 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  Vector<Element> elements_regular;
  Vector<Element> elements_cohesive;

  splitElementByKind(elements, elements_regular, elements_cohesive);

  SolidMechanicsModel::packElementData(buffer, elements_regular, tag);

#ifndef AKANTU_NDEBUG
  packBarycenter(buffer, elements_cohesive);
#endif

  switch(tag) {

  case _gst_material_id: {
    packElementalDataHelper(element_index_by_material, buffer, elements_cohesive, false, "CohesiveFEM");
    break;
  }
  case _gst_smm_boundary: {
    packNodalDataHelper(*force, buffer, elements_cohesive);
    packNodalDataHelper(*velocity, buffer, elements_cohesive);
    packNodalDataHelper(*boundary, buffer, elements_cohesive);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Vector<Element> * elements_per_mat = new Vector<Element>[materials.size()];
    splitElementByMaterial(elements_cohesive, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->packElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelCohesive::unpackElementData(CommunicationBuffer & buffer,
							   const Vector<Element> & elements,
							   SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  Vector<Element> elements_regular;
  Vector<Element> elements_cohesive;

  splitElementByKind(elements, elements_regular, elements_cohesive);

  SolidMechanicsModel::unpackElementData(buffer, elements_regular, tag);

#ifndef AKANTU_NDEBUG
  unpackBarycenter(buffer, elements_cohesive, tag);
#endif

  switch(tag) {
  case _gst_material_id: {
    unpackElementalDataHelper(element_index_by_material, buffer, elements_cohesive, false, "CohesiveFEM");
    break;
  }
  case _gst_smm_boundary: {
    unpackNodalDataHelper(*force, buffer, elements_cohesive);
    unpackNodalDataHelper(*velocity, buffer, elements_cohesive);
    unpackNodalDataHelper(*boundary, buffer, elements_cohesive);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Vector<Element> * elements_per_mat = new Vector<Element>[materials.size()];
    splitElementByMaterial(elements_cohesive, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->unpackElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
}
