/**
 * @file   historic_velocity_fric_coef.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Nov  1 14:34:01 2011
 *
 * @brief  implementation of historic velocity dependent friction coefficient
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
#include "historic_velocity_fric_coef.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
HistoricVelocityFricCoef::HistoricVelocityFricCoef(ContactRigid & contact,
						   const Surface & master_surface,
						   const Real beta) :
  FrictionCoefficient(contact, master_surface), 
  spatial_dimension(contact.getSpatialDimension()), 
  beta(beta),  
  nodes(0,1),
  active(0,1),
  contact_time(0,1) {
  AKANTU_DEBUG_IN();

  // compute time range of history
  Real tolerance_epsilon = 1e-5;
  Real total_time = log(1/tolerance_epsilon) / beta;
  Real time_step = this->contact.getModel().getTimeStep();
  UInt nb_time_steps = (UInt)ceil(total_time / time_step) + 1;

  std::cout << " * Historic Velocity Friction Coefficient: total history time = " << total_time << " with numbers of time steps = " << nb_time_steps << std::endl;

  // declare vector
  this->weights = new Vector<Real>(nb_time_steps,1);
  this->historic_velocities = new CircularVector< Vector<Real> * >(nb_time_steps, 1);

  // fill the weights vector
  Real time = 0.;
  for (Int i=nb_time_steps-1; i>=0; --i) {
    (*weights)(i) = exp(-beta*time);
    time += time_step;
  }

  this->generalized_sliding_velocities = new Vector<Real>(0,1);
  this->node_stick_status = new Vector<bool>(0,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
HistoricVelocityFricCoef::HistoricVelocityFricCoef(ContactRigid & contact,
						   const Surface & master_surface) :
  FrictionCoefficient(contact, master_surface), 
  spatial_dimension(contact.getSpatialDimension()), 
  beta(1.),  
  nodes(0,1),
  active(0,1),
  contact_time(0,1) {
  AKANTU_DEBUG_IN();

  std::cout << " * Historic Velocity Friction Coefficient: not used, only inherited from." << std::endl;

  // declare vector
  this->weights = new Vector<Real>(0,1);
  this->historic_velocities = new CircularVector< Vector<Real> * >(0, 1);
  this->generalized_sliding_velocities = new Vector<Real>(0,1);
  this->node_stick_status = new Vector<bool>(0,1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
HistoricVelocityFricCoef::~HistoricVelocityFricCoef() {
  AKANTU_DEBUG_IN();
  
  delete this->weights;
  
  for (UInt i=0; i<historic_velocities->getSize(); ++i)
    delete ((*historic_velocities)(i));
  delete this->historic_velocities;

  delete this->generalized_sliding_velocities;
  delete this->node_stick_status;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HistoricVelocityFricCoef::initializeComputeFricCoef() {
  AKANTU_DEBUG_IN();
  
  /*
    Find if more nodes have to be considered
   */
  
  // get access to impactors information 
  const ContactRigid::SurfaceToImpactInfoMap & imp_info = this->contact.getImpactorsInformation();
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it;
  
  it = imp_info.find(this->master_surface);
  AKANTU_DEBUG_ASSERT(it != imp_info.end(), 
		      "Couldn't find impactor information object for master surface " << master_surface);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it->second;
  Vector<UInt> * active_nodes = impactor_info->active_impactor_nodes;
  UInt nb_active_nodes = active_nodes->getSize();

  Real time_step = this->contact.getModel().getTimeStep();

  // check if in nodes vector. if no, put it.
  this->active.clear();
  for (UInt i=0; i<nb_active_nodes; ++i) {
    UInt node = (*active_nodes)(i);
    Int pos = this->nodes.find(node);
    if (pos == -1) {
      nodes.push_back(node);
      active.push_back(true);
      contact_time.push_back(0.);
    }
    else {
      active(pos) = true;
      contact_time(pos) += time_step;
    }
  }
  UInt nb_nodes = this->nodes.getSize();
  
  // resize the vectors in the historic velocities
  UInt nb_historic_velocities = this->historic_velocities->getSize();
  for (UInt i=0; i<nb_historic_velocities; ++i) {
    Vector<Real> * historic_velocity = (*historic_velocities)(i);
    if (historic_velocity != NULL) {
      historic_velocity->resize(nb_nodes); // the new values are initialized to zero
    }
  }

  /* 
     find the velocities of this time step and put in vector
   */

  // advance circular vector and access the velocities
  historic_velocities->makeStep();
  if ((*historic_velocities)(nb_historic_velocities-1) == NULL)
     (*historic_velocities)(nb_historic_velocities-1) =  new Vector<Real>(nb_nodes,1);
  else {
    (*historic_velocities)(nb_historic_velocities-1)->clear();
  }
  Vector<Real> * current_velocities = (*historic_velocities)(nb_historic_velocities-1);
  Real * velocity_val = this->contact.getModel().getVelocity().values;
  
  // access the master normals
  Vector<Real> * master_normals = impactor_info->master_normals;
  Real * master_normals_val = master_normals->values;

  // norms for each point
  Vector<Real> norms(nb_nodes,1);

  // loop over all nodes in here
  for(UInt n = 0; n < nb_nodes; ++n) {
    // if node is not active set velocity to zero
    if (!active(n)) {
      (*current_velocities)(n) = 0.;
      contact_time(n) = 0.; // put contact_time to zero for non active nodes
    }
    
    // OTHERWISE ...
    else {

      // ... COMPUTE TANGENTIAL
      UInt impactor = nodes(n);
      
      // get indexes in the vectors
      UInt master_normal_index = active_nodes->find(impactor) * this->spatial_dimension;
      UInt velocity_index = impactor * this->spatial_dimension;
      
      // project slave velocity
      Real projected_length;
      if(this->spatial_dimension == 2)
	projected_length = Math::vectorDot2(&master_normals_val[master_normal_index], &velocity_val[velocity_index]);
      else if(this->spatial_dimension == 3)
	projected_length = Math::vectorDot3(&master_normals_val[master_normal_index], &velocity_val[velocity_index]);
      
      Real tangential_impactor_velocity[3];
      for(UInt i=0; i<this->spatial_dimension; ++i) {
	tangential_impactor_velocity[i] = velocity_val[velocity_index + i] - projected_length * master_normals_val[master_normal_index + i];
      }
      
      // get master velocity
      Real tangential_master_velocity[3];
      computeTangentialMasterVelocity(n, impactor_info, tangential_master_velocity);
      
      // compute relative sliding speed
      Real relative_sliding_speed;
      Real relative_sliding_velocity[3];
      for(UInt i =0; i<this->spatial_dimension; ++i)
	relative_sliding_velocity[i] = tangential_impactor_velocity[i] - tangential_master_velocity[i];
      if(this->spatial_dimension == 2)
	relative_sliding_speed = Math::norm2(relative_sliding_velocity);
      else if(this->spatial_dimension == 3)
	relative_sliding_speed = Math::norm3(relative_sliding_velocity);
      
      (*current_velocities)(n) = relative_sliding_speed;

      // ... COMPUTE NORM
      if (contact_time(n)/time_step < 1e-14)
	norms(n) = 1.; // there should be only the instant velocity
      else {
	norms(n) = beta / (1 - exp(-beta*contact_time(n)));
      }
    }
  }
  
  /*
    Compute the weighted average of sliding velocity
  */
  
  Vector<Real> results(nb_nodes, 1);
  for (UInt i=0; i<nb_historic_velocities; ++i) {
    Vector<Real> * velocities = (*historic_velocities)(i);
    if (velocities == NULL) continue;

    // set inactive nodes to zero
    for (UInt j=0; j<nb_nodes; ++j) {
      if (!active(j))
	(*velocities)(j) = 0.;
    }

    // multiply with weight   
    Real weight = (*weights)(i);
    Vector<Real> weighted_velocities(0,1);
    weighted_velocities.copy(*velocities);
    weighted_velocities *= weight * time_step;
    
    // add up
    results += weighted_velocities;
  }

  // divide by normalized velocity
  for (UInt j=0; j<nb_nodes; ++j) {
    results(j) *= norms(j);
  }

  /*
    fill the final vector
  */
  
  // resize the relative sliding velocity vector to the nb of active impactor nodes
  this->generalized_sliding_velocities->resize(nb_active_nodes);
  this->node_stick_status->resize(nb_active_nodes);
  Vector<bool> * stick_info = impactor_info->node_is_sticking;

  for (UInt i=0; i<nb_active_nodes; ++i) {
    UInt node = (*active_nodes)(i);
    UInt pos = nodes.find(node);
    (*generalized_sliding_velocities)(i) = results(pos);
    (*node_stick_status)(i) = (*stick_info)(i,0);
    //  std::cout << (*generalized_sliding_velocities)(i) << " ";
  }
  //std::cout << std::endl;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void HistoricVelocityFricCoef::computeTangentialMasterVelocity(UInt impactor_index, 
							       ContactRigid::ImpactorInformationPerMaster * impactor_info, 
							       Real * tangential_master_velocity) {
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
    case _quadrangle_8:
    case _hexahedron_8:
    case _point: 
    case _bernoulli_beam_2:
    case _max_element_type: {
      AKANTU_DEBUG_ERROR("Not a valid surface element type : " << type << " for computation of tangential velocity of master element");
      break;
    }
    }
  }
  
  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
