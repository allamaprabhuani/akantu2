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

    UInt mat = element_material(element.type, _not_ghost)(element.element);
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  default: {
    UInt mat = element_material(element.type, _ghost)(element.element);
    size += materials[mat]->getNbDataToPack(element, tag);
    //AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
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

    UInt mat = element_material(element.type, _ghost)(element.element);
    size += materials[mat]->getNbDataToPack(element, tag);
    break;
  }
  case _gst_smm_boundary: {
    // force, displacement, boundary
    size += nb_nodes_per_element * spatial_dimension * (2 * sizeof(Real) + sizeof(bool));
    break;
  }
  default: {
    UInt mat = element_material(element.type, _ghost)(element.element);
    size += materials[mat]->getNbDataToPack(element, tag);
    // AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
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
  UInt * conn  = mesh.getConnectivity(element.type, ghost_type).storage();

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

    UInt mat = element_material(element.type, ghost_type)(element.element);
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
    UInt mat = element_material(element.type, ghost_type)(element.element);
    materials[mat]->packData(buffer, element, tag);
    //AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModel::unpackData(CommunicationBuffer & buffer,
					    const Element & element,
					    SynchronizationTag tag) {
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

    UInt mat = element_material(element.type, ghost_type)(element.element);
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
    UInt mat = element_material(element.type, ghost_type)(element.element);
    materials[mat]->unpackData(buffer, element, tag);
    //AKANTU_DEBUG_ERROR("Unknown ghost synchronization tag : " << tag);
  }
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
    Vector<Real>::iterator<types::RVector> it_disp = displacement->begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> it_velo = velocity->begin(spatial_dimension);
    buffer << it_disp[index];
    buffer << it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Vector<Real>::iterator<types::RVector> it_res = residual->begin(spatial_dimension);
    buffer << it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("pack mass of node " << index << " which is " << (*mass)(index,0));
    Vector<Real>::iterator<types::RVector> it_mass = mass->begin(spatial_dimension);
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
    Vector<Real>::iterator<types::RVector> it_disp = displacement->begin(spatial_dimension);
    Vector<Real>::iterator<types::RVector> it_velo = velocity->begin(spatial_dimension);
    buffer >> it_disp[index];
    buffer >> it_velo[index];
    break;
  }
  case _gst_smm_res: {
    Vector<Real>::iterator<types::RVector> it_res = residual->begin(spatial_dimension);
    buffer >> it_res[index];
    break;
  }
  case _gst_smm_mass: {
    AKANTU_DEBUG_INFO("mass of node " << index << " was " << (*mass)(index,0));
    Vector<Real>::iterator<types::RVector> it_mass = mass->begin(spatial_dimension);
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
void SolidMechanicsModel::solveDynamic(Vector<Real> & increment) {
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
