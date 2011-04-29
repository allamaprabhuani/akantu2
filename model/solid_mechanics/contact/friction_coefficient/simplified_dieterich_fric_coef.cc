/**
 * @file   simplified_dieterich_fric_coef.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Fri Mar  4 15:00:32 2011
 *
 * @brief  implementation of the simpliefied dieterich friction law
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
#include "simplified_dieterich_fric_coef.hh"
#include "contact_rigid.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
SimplifiedDieterichFricCoef::SimplifiedDieterichFricCoef(ContactRigid & contact,
							 const Surface & master_surface) : FrictionCoefficient(contact, master_surface), spatial_dimension(this->contact.getSpatialDimension()), mu_zero(0.), a_factor(0.), b_factor(0.), v_normalizer(1.), theta_normalizer(1.) {
  AKANTU_DEBUG_IN();

  this->relative_sliding_velocities = new Vector<Real>(0,1);
  this->theta_state_variables = new Vector<Real>(0,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
SimplifiedDieterichFricCoef::~SimplifiedDieterichFricCoef() {
  AKANTU_DEBUG_IN();
  
  delete this->relative_sliding_velocities;
  delete this->theta_state_variables;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SimplifiedDieterichFricCoef::initializeComputeFricCoef() {
  AKANTU_DEBUG_IN();

  // find impactor_information for given master
  const ContactRigid::SurfaceToImpactInfoMap & imp_info = this->contact.getImpactorsInformation();
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it;

  it = imp_info.find(this->master_surface);
  AKANTU_DEBUG_ASSERT(it != imp_info.end(), 
		      "Couldn't find impactor information object for master surface " << master_surface);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it->second;

  /*
  const std::vector<ContactRigid::ImpactorInformationPerMaster * > impactor_information_vector = this->contact.getImpactorsInformation();

  // find index of master surface in impactors_information 
  Int master_index = -1;
  for (UInt i=0; i < impactor_information_vector.size(); ++i) {
    if (impactor_information_vector.at(i)->master_id == this->master_surface) {
      master_index = i;
      break;
    }
  }
  AKANTU_DEBUG_ASSERT(master_index != -1, "No impactor information object for master" << this->master_surface << "in impactors_information vector that is ");

  ContactRigid::ImpactorInformationPerMaster * impactor_info = impactor_information_vector.at(master_index);
  */

  Vector<UInt> * active_nodes = impactor_info->active_impactor_nodes;
  UInt * active_nodes_val = active_nodes->values;
  Vector<Real> * master_normals = impactor_info->master_normals;
  Real * master_normals_val = master_normals->values;

  Real * velocity_val = this->contact.getModel().getVelocity().values;

  // resize the relative sliding velocity vector to the nb of active impactor nodes
  this->relative_sliding_velocities->resize(active_nodes->getSize());
  Real * relative_sliding_velocities_val = this->relative_sliding_velocities->values;

  for(UInt n = 0; n < active_nodes->getSize(); ++n) {
    UInt impactor = active_nodes_val[n];
    UInt master_normal_index = n * this->spatial_dimension;
    UInt velocity_index = impactor * this->spatial_dimension;
    
    Real projected_length;
    if(this->spatial_dimension == 2)
      projected_length = Math::vectorDot2(&master_normals_val[master_normal_index], &velocity_val[velocity_index]);
    else if(this->spatial_dimension == 3)
      projected_length = Math::vectorDot3(&master_normals_val[master_normal_index], &velocity_val[velocity_index]);

    Real tangential_impactor_velocity[3];
    for(UInt i=0; i<this->spatial_dimension; ++i) {
      tangential_impactor_velocity[i] = velocity_val[velocity_index + i] - projected_length * master_normals_val[master_normal_index + i];
    }

    Real tangential_master_velocity[3];
    //Real * tangential_master_velocity_p = &(tangential_master_velocity[0]);
    computeTangentialMasterVelocity(n, impactor_info, tangential_master_velocity);
    
    Real relative_sliding_speed;
    Real relative_sliding_velocity[3];
    //Real * relative_sliding_velocity_p = &(relative_sliding_velocity_p[0]);
    for(UInt i =0; i<this->spatial_dimension; ++i)
      relative_sliding_velocity[i] = tangential_impactor_velocity[i] - tangential_master_velocity[i];
    if(this->spatial_dimension == 2)
      relative_sliding_speed = Math::norm2(relative_sliding_velocity);
    else if(this->spatial_dimension == 3)
      relative_sliding_speed = Math::norm3(relative_sliding_velocity);
    
    relative_sliding_velocities_val[n] = relative_sliding_speed; 

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SimplifiedDieterichFricCoef::computeTangentialMasterVelocity(UInt impactor_index, ContactRigid::ImpactorInformationPerMaster * impactor_info, Real * tangential_master_velocity) {
  AKANTU_DEBUG_IN();

  // for ContactRigid, the master velocity is by definition zero
  const ContactType contact_type = this->contact.getType();
  if(contact_type == _ct_rigid /* do not change for general contact */) {
    for (UInt i=0; i<this->spatial_dimension; ++i)
      tangential_master_velocity[i] = 0.;
  }
  else {

    // compute the tangential master velocity
    ElementType type = impactor_info->master_element_type->at(impactor_index);
    switch(type) {
      /*
	case _segment_2: {
	//computeComponentsOfProjectionSegment2(impactor_node, surface_element, normal, gap, projected_position);
	break;
	}
      */
      /*
	case _triangle_3: {
	//computeComponentsOfProjectionTriangle3(impactor_node, surface_element, normal, gap, projected_position);
	break;
	}
      */
    case _not_defined: {
      AKANTU_DEBUG_ERROR("Not a valid surface element type : " << type << " for computation of tangential velocity of master element");
      break;
    }
    case _segment_2:
    case _segment_3:
    case _triangle_3:
    case _triangle_6:
    case _tetrahedron_4:
    case _tetrahedron_10:
    case _quadrangle_4:
    case _hexahedron_8:
    case _point: 
    case _max_element_type: {
      AKANTU_DEBUG_ERROR("Not a valid surface element type : " << type << " for computation of tangential velocity of master element");
      break;
    }
    }
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SimplifiedDieterichFricCoef::setParam(const std::string & key, const std::string & value) {
  AKANTU_DEBUG_IN();
  
  std::stringstream sstr(value);
  if(key == "mu_zero") { sstr >> this->mu_zero; }
  else if(key == "a_factor") { sstr >> this->a_factor; }
  else if(key == "b_factor") { sstr >> this->b_factor; }
  else if(key == "v_normalizer") { sstr >> this->v_normalizer; }
  else if(key == "theta_normalizer") { sstr >> this->theta_normalizer; }
  else { FrictionCoefficient::setParam(key, value); }
  
  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
