/**
 * @file   solid_mechanics_model_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 04 10:58:42 2010
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
inline FEM & SolidMechanicsModel::getFEMBoundary(const ID & name) {
  return dynamic_cast<FEM &>(getFEMClassBoundary<MyFEMType>(name));
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::splitElementByMaterial(const Array<Element> & elements,
						 Array<Element> * elements_per_mat) const {
  ElementType current_element_type = _not_defined;
  GhostType current_ghost_type = _casper;
  UInt * elem_mat = NULL;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    if(el.type != current_element_type || el.ghost_type != current_ghost_type) {
      current_element_type = el.type;
      current_ghost_type   = el.ghost_type;
      elem_mat = element_index_by_material(el.type, el.ghost_type).storage();
    }
    elements_per_mat[elem_mat[2*el.element+1]].push_back(el);
  }
}

/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataForElements(const Array<Element> & elements,
						      SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;

#ifndef AKANTU_NDEBUG
  size += elements.getSize() * spatial_dimension * sizeof(Real); /// position of the barycenter of the element (only for check)
#endif

  UInt nb_nodes_per_element = 0;

  Array<Element>::const_iterator<Element> it  = elements.begin();
  Array<Element>::const_iterator<Element> end = elements.end();
  for (; it != end; ++it) {
    const Element & el = *it;
    nb_nodes_per_element += Mesh::getNbNodesPerElement(el.type);
  }

  switch(tag) {
  case _gst_material_id: {
    size += elements.getSize() * 2 * sizeof(UInt);
    break;
  }
  case _gst_smm_mass: {
    size += nb_nodes_per_element * sizeof(Real) * spatial_dimension; // mass vector
    break;
  }
  case _gst_smm_for_strain: {
    size += nb_nodes_per_element * spatial_dimension * sizeof(Real); // displacement
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
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    this->splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      size += materials[i]->getNbDataForElements(elements_per_mat[i], tag);
    }
    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::packBarycenter(CommunicationBuffer & buffer,
						const Array<Element> & elements,
                                                SynchronizationTag tag) const {
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();
  for (; bit != bend; ++bit) {
    const Element & element = *bit;
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter.storage(), element.ghost_type);
    buffer << barycenter;
  }
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::packElementData(CommunicationBuffer & buffer,
						 const Array<Element> & elements,
						 SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  packBarycenter(buffer, elements, tag);
#endif

  switch(tag) {

  case _gst_material_id: {
    packElementalDataHelper(element_index_by_material, buffer, elements, false);
    break;
  }
  case _gst_smm_mass: {
    packNodalDataHelper(*mass, buffer, elements);
    break;
  }
  case _gst_smm_for_strain: {
    packNodalDataHelper(*displacement, buffer, elements);
    break;
  }
  case _gst_smm_boundary: {
    packNodalDataHelper(*force, buffer, elements);
    packNodalDataHelper(*velocity, buffer, elements);
    packNodalDataHelper(*boundary, buffer, elements);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->packElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackBarycenter(CommunicationBuffer & buffer,
						  const Array<Element> & elements,
						  SynchronizationTag tag) {
  Array<Element>::const_iterator<Element> bit  = elements.begin();
  Array<Element>::const_iterator<Element> bend = elements.end();

  for (; bit != bend; ++bit) {
    const Element & element = *bit;

    Vector<Real> barycenter_loc(spatial_dimension);
    mesh.getBarycenter(element.element, element.type, barycenter_loc.storage(), element.ghost_type);
    Vector<Real> barycenter(spatial_dimension);
    buffer >> barycenter;
    Real tolerance = 1e-15;
    for (UInt i = 0; i < spatial_dimension; ++i) {
      if(!(std::abs((barycenter(i) - barycenter_loc(i))/barycenter_loc(i)) <= tolerance))
	AKANTU_DEBUG_ERROR("Unpacking an unknown value for the element: "
			   << element
			   << "(barycenter[" << i << "] = " << barycenter_loc(i)
			   << " and buffer[" << i << "] = " << barycenter(i) << ") ["
			   << std::abs((barycenter(i) - barycenter_loc(i))/barycenter_loc(i))
			   << "] - tag: " << tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackElementData(CommunicationBuffer & buffer,
						   const Array<Element> & elements,
						   SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  unpackBarycenter(buffer, elements, tag);
#endif

  switch(tag) {
  case _gst_material_id: {
    unpackElementalDataHelper(element_index_by_material, buffer, elements, false);
    break;
  }
  case _gst_smm_mass: {
    unpackNodalDataHelper(*mass, buffer, elements);
    break;
  }
  case _gst_smm_for_strain: {
    unpackNodalDataHelper(*displacement, buffer, elements);
    break;
  }
  case _gst_smm_boundary: {
    unpackNodalDataHelper(*force, buffer, elements);
    unpackNodalDataHelper(*velocity, buffer, elements);
    unpackNodalDataHelper(*boundary, buffer, elements);
    break;
  }
  default: {
  }
  }

  if(tag != _gst_material_id) {
    Array<Element> * elements_per_mat = new Array<Element>[materials.size()];
    splitElementByMaterial(elements, elements_per_mat);

    for (UInt i = 0; i < materials.size(); ++i) {
      materials[i]->unpackElementData(buffer, elements_per_mat[i], tag);
    }

    delete [] elements_per_mat;
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
inline UInt SolidMechanicsModel::getNbDataToPack(SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  //  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_res: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_smm_mass: {
    size += sizeof(Real) * spatial_dimension;
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
  //  UInt nb_nodes = mesh.getNbNodes();

  switch(tag) {
  case _gst_smm_uv: {
    size += sizeof(Real) * spatial_dimension * 2;
    break;
  }
  case _gst_smm_res: {
    size += sizeof(Real) * spatial_dimension;
    break;
  }
  case _gst_smm_mass: {
    size += sizeof(Real) * spatial_dimension;
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
    Array<Real>::iterator< Vector<Real> > it_disp = displacement->begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > it_velo = velocity->begin(spatial_dimension);
    buffer << it_disp[index];
    buffer << it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Array<Real>::iterator< Vector<Real> > it_res = residual->begin(spatial_dimension);
    buffer << it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("pack mass of node " << index << " which is " << (*mass)(index,0));
    Array<Real>::iterator< Vector<Real> > it_mass = mass->begin(spatial_dimension);
    buffer << it_mass[index];
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
					    SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  switch(tag) {
  case _gst_smm_uv: {
    Array<Real>::iterator< Vector<Real> > it_disp = displacement->begin(spatial_dimension);
    Array<Real>::iterator< Vector<Real> > it_velo = velocity->begin(spatial_dimension);
    buffer >> it_disp[index];
    buffer >> it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Array<Real>::iterator< Vector<Real> > it_res = residual->begin(spatial_dimension);
    buffer >> it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("mass of node " << index << " was " << (*mass)(index,0));
    Array<Real>::iterator< Vector<Real> > it_mass = mass->begin(spatial_dimension);
    buffer >> it_mass[index];
    AKANTU_DEBUG_INFO("mass of node " << index << " is now " << (*mass)(index,0));
    break;
  }
  default: {
    AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
#include "sparse_matrix.hh"
#include "solver.hh"
__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<NewmarkBeta::IntegrationSchemeCorrectorType type>
void SolidMechanicsModel::solveDynamic(Array<Real> & increment) {
  AKANTU_DEBUG_INFO("Solving Ma + Cv + Ku = f");

  NewmarkBeta * nmb_int = dynamic_cast<NewmarkBeta *>(integrator);
  Real c = nmb_int->getAccelerationCoefficient<type>(time_step);
  Real d = nmb_int->getVelocityCoefficient<type>(time_step);
  Real e = nmb_int->getDisplacementCoefficient<type>(time_step);

  // A = c M + d C + e K
  jacobian_matrix->clear();

  if(type != NewmarkBeta::_acceleration_corrector)
    jacobian_matrix->add(*stiffness_matrix, e);

  jacobian_matrix->add(*mass_matrix, c);

  mass_matrix->saveMatrix("M.mtx");
  if(velocity_damping_matrix)
    jacobian_matrix->add(*velocity_damping_matrix, d);

  jacobian_matrix->applyBoundary(*boundary);

#ifndef AKANTU_NDEBUG
  if(AKANTU_DEBUG_TEST(dblDump))
    jacobian_matrix->saveMatrix("J.mtx");
#endif

  jacobian_matrix->saveMatrix("J.mtx");

  solver->setRHS(*residual);

  // solve A w = f
  solver->solve(increment);
}
/* -------------------------------------------------------------------------- */
