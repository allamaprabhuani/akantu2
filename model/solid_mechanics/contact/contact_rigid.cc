/**
 * @file   contact_rigid.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 17:15:08 2010
 *
 * @brief Specialization  of the  contact structure for  3d contact  in explicit
 * time scheme
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
#include "contact_rigid.hh"
#include "contact_search.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ContactRigid::ImpactorInformationPerMaster::ImpactorInformationPerMaster(const Surface master_id, const UInt spatial_dimension) : master_id(master_id) {
  AKANTU_DEBUG_IN();

  this->impactor_surfaces = new std::vector<Surface>(0,1);

  this->master_normals = new Vector<Int>(0, spatial_dimension);
  this->active_impactor_nodes = new Vector<UInt>(0,1);
  this->node_is_sticking = new Vector<bool>(0,2);
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactRigid::ImpactorInformationPerMaster::~ImpactorInformationPerMaster() {
  AKANTU_DEBUG_IN();

  delete this->impactor_surfaces;  

  delete this->master_normals;
  delete this->active_impactor_nodes;
  delete this->node_is_sticking;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactRigid::ContactRigid(const SolidMechanicsModel & model,
			   const ContactType & type,
			   const ContactID & id,
			   const MemoryID & memory_id) :
  Contact(model, type, id, memory_id), spatial_dimension(model.getSpatialDimension()), mesh(model.getFEM().getMesh()) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactRigid::~ContactRigid() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  Contact::addMasterSurface(master_surface);

  ImpactorInformationPerMaster * impactor_info = new ImpactorInformationPerMaster(master_surface, this->spatial_dimension);
  impactors_information.push_back(impactor_info);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addImpactorSurfaceToMasterSurface(const Surface & impactor_surface, const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  
  for (UInt m=0; m < this->impactors_information.size(); ++m) {
    if (this->impactors_information.at(m)->master_id == master_surface) {
      this->impactors_information.at(m)->impactor_surfaces->push_back(impactor_surface);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::removeMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  Contact::removeMasterSurface(master_surface);

  for (UInt m=0; m < this->impactors_information.size(); ++m) {
    if (this->impactors_information.at(m)->master_id == master_surface) {
      this->impactors_information.erase(this->impactors_information.begin()+m);
      break;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::removeImpactorSurfaceFromMasterSurface(const Surface & impactor_surface, const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  
  for (UInt m=0; m < this->impactors_information.size(); ++m) {
    ImpactorInformationPerMaster * impactor_info = this->impactors_information.at(m);
    if (impactor_info->master_id == master_surface) {
      for (UInt i=0; i < impactor_info->impactor_surfaces->size(); ++i) {
	Surface imp_surface = impactor_info->impactor_surfaces->at(i);
	if (imp_surface == impactor_surface) {
	  impactor_info->impactor_surfaces->erase(impactor_info->impactor_surfaces->begin()+i);
	  break;
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::solveContact() {
  AKANTU_DEBUG_IN();

  for(UInt m=0; m < master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);
    PenetrationList * penet_list = new PenetrationList();
    contact_search->findPenetration(master, *penet_list);
    
    // find index of master surface in impactors_information 
    Int master_index = -1;
    for (UInt i=0; i < this->impactors_information.size(); ++i) {
      if (this->impactors_information.at(i)->master_id == master) {
	master_index = i;
	break;
      }
    }
    AKANTU_DEBUG_ASSERT(master_index != -1, "No impactor information object for master" << master << "in impactors_information vector that is ");
    
    solvePenetrationClosestProjection(master_index, *penet_list);
    delete penet_list;
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*void ContactRigid::solvePenetration(const PenetrationList & penet_list) {
  AKANTU_DEBUG_IN();

  const UInt dim = ;

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  /// find existing surface element types
  UInt nb_types = type_list.size();
  UInt nb_facet_types = 0;
  ElementType facet_type[_max_element_type];
  
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(mesh.getSpatialDimension(type) == spatial_dimension) {
      ElementType current_facet_type = mesh.getFacetElementType(type);
      facet_type[nb_facet_types++] = current_facet_type;
    }
  }

  Real * current_position = model.getCurrentPosition().values;
  Real * displacement     = model.getDisplacement().values;
  Real * increment        = model.getIncrement().values;

  UInt * penetrating_nodes = penet_list.penetrating_nodes.values;

  for (UInt n=0; n < penet_list.penetrating_nodes.getSize(); ++n) {
    UInt current_node = penetrating_nodes[n];
        
    // count number of elements penetrated by this node
    UInt nb_penetrated_facets = 0;
    ElementType penetrated_type;
    for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {
      ElementType type = facet_type[el_type];
      UInt offset_min = penet_list.penetrated_facets_offset[type].get(n);
      UInt offset_max = penet_list.penetrated_facets_offset[type].get(n+1);
      nb_penetrated_facets += offset_max - offset_min;
      penetrated_type = type;
    }

    // easy case: node penetrated only one facet
    if(nb_penetrated_facets == 1) {
      Real * facets_normals      = penet_list.facets_normals[penetrated_type].values;
      Real * gaps                = penet_list.gaps[penetrated_type].values;
      Real * projected_positions = penet_list.projected_positions[penetrated_type].values;

      UInt offset_min = penet_list.penetrated_facets_offset[penetrated_type].get(n);
      for(UInt i=0; i < dim; ++i) {
	current_position[current_node*dim + i] = projected_positions[offset_min*dim + i];
	Real displacement_correction = gaps[offset_min*dim + i] * normals[offset_min*dim + i];
	displacement[current_node*dim + i] = displacement[current_node*dim + i] - displacement_correction;
	increment   [current_node*dim + i] = increment   [current_node*dim + i] - displacement_correction;
      }
    }
    
    // if more penetrated facets, find the one from which the impactor node is coming
    else {

    }
  }
  }*/

/* -------------------------------------------------------------------------- */
void ContactRigid::solvePenetrationClosestProjection(const UInt master_index, 
						     const PenetrationList & penet_list) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  /// find existing surface element types
  UInt nb_types = type_list.size();
  UInt nb_facet_types = 0;
  ElementType facet_type[_max_element_type];
  
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    ElementType type = *it;
    if(mesh.getSpatialDimension(type) == spatial_dimension) {
      ElementType current_facet_type = mesh.getFacetElementType(type);
      facet_type[nb_facet_types++] = current_facet_type;
    }
  }

  for (UInt n=0; n < penet_list.penetrating_nodes.getSize(); ++n) {
            
    // find facet for which the gap is the smallest
    Real min_gap = std::numeric_limits<Real>::max();
    ElementType penetrated_type;
    UInt penetrated_facet_offset;
    for (UInt el_type = 0; el_type < nb_facet_types; ++el_type) {
      ElementType type = facet_type[el_type];
      Real * gaps = penet_list.gaps[type]->values;
      UInt offset_min = penet_list.penetrated_facets_offset[type]->get(n);
      UInt offset_max = penet_list.penetrated_facets_offset[type]->get(n+1);
      for (UInt f = offset_min; f < offset_max; ++f) {
	if(gaps[f] < min_gap) {
	  min_gap = gaps[f];
	  penetrated_type = type;
	  penetrated_facet_offset = f;
	}
      }
    }

    // correct the position of the impactor
    projectImpactor(penet_list, n, penetrated_type, penetrated_facet_offset);
    lockImpactorNode(master_index, penet_list, n, penetrated_type, penetrated_facet_offset);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::projectImpactor(const PenetrationList & penet_list, const UInt impactor_index, const ElementType facet_type, const UInt facet_offset) {

  AKANTU_DEBUG_IN();
  
  const bool increment_flag = model.getIncrementFlag();

  UInt * penetrating_nodes   = penet_list.penetrating_nodes.values;
  Real * facets_normals      = penet_list.facets_normals[facet_type]->values;
  Real * gaps                = penet_list.gaps[facet_type]->values;
  Real * projected_positions = penet_list.projected_positions[facet_type]->values;

  Real * current_position = model.getCurrentPosition().values;
  Real * displacement     = model.getDisplacement().values;
  Real * increment = NULL;
  if(increment_flag)
    increment = model.getIncrement().values;

  UInt impactor_node = penetrating_nodes[impactor_index];

  for(UInt i=0; i < spatial_dimension; ++i) {
    current_position[impactor_node*spatial_dimension + i] = projected_positions[facet_offset*spatial_dimension + i];
    Real displacement_correction = gaps[facet_offset] * facets_normals[facet_offset*spatial_dimension + i];
    displacement[impactor_node*spatial_dimension + i] = displacement[impactor_node*spatial_dimension + i] - displacement_correction;
    if(increment_flag)
      increment   [impactor_node*spatial_dimension + i] = increment   [impactor_node*spatial_dimension + i] - displacement_correction;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::lockImpactorNode(const UInt master_index, const PenetrationList & penet_list, const UInt impactor_index, const ElementType facet_type, const UInt facet_offset) {
  AKANTU_DEBUG_IN();
  
  UInt * penetrating_nodes = penet_list.penetrating_nodes.values;
  UInt impactor_node = penetrating_nodes[impactor_index];

  Real * facets_normals = penet_list.facets_normals[facet_type]->values;
  Real * facet_normal = &facets_normals[facet_offset*spatial_dimension];
  Int normal[this->spatial_dimension];
  //Int * normal_val = &normal[0];

  bool * bound_val = this->model.getBoundary().values;
  Real * veloc_val = this->model.getVelocity().values;
  Real * accel_val = this->model.getAcceleration().values;

  for(UInt i = 0; i < spatial_dimension; ++i)
    normal[i] = static_cast<Int>(floor(facet_normal[i] + 0.5));
  
  for(UInt i = 0; i < spatial_dimension; ++i) {
    if(normal[i] != 0) {
      UInt index = impactor_node * spatial_dimension + i;
      bound_val[index] = true;
      veloc_val[index] = 0.;
      accel_val[index] = 0.;
    }
  }

  ImpactorInformationPerMaster * impactor_info = this->impactors_information.at(master_index);
  impactor_info->active_impactor_nodes->push_back(impactor_node);
  impactor_info->master_normals->push_back(normal);
  Real init_sticking[2];
  init_sticking[0] = true; init_sticking[1] = true;
  impactor_info->node_is_sticking->push_back(init_sticking);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::avoidAdhesion() {
  AKANTU_DEBUG_IN();

  Real * residual_val = this->model.getResidual().values;
  bool * bound_val    = this->model.getBoundary().values;

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    // find index of master surface in impactors_information 
    Int master_index = -1;
    for (UInt i=0; i < this->impactors_information.size(); ++i) {
      if (this->impactors_information.at(i)->master_id == master) {
	master_index = i;
	break;
      }
    }
    AKANTU_DEBUG_ASSERT(master_index != -1, "No impactor information object for master" << master << "in impactors_information vector that is ");

    ImpactorInformationPerMaster * impactor_info = this->impactors_information.at(master_index);
    
    for (UInt n=0; n < impactor_info->active_impactor_nodes->getSize(); ++n) {
      UInt current_node = impactor_info->active_impactor_nodes->at(n);
      
      for (UInt i=0; i < spatial_dimension; ++i) {
	Int direction = impactor_info->master_normals->at(n,i);
	Real force = residual_val[current_node * spatial_dimension + i];
	if(force * direction > 0.) {
	  bound_val[current_node * spatial_dimension + i] = false;
	  impactor_info->active_impactor_nodes->erase(n);
	  impactor_info->master_normals->erase(n);
	  impactor_info->node_is_sticking->erase(n);
	  n--;
	  break;
	}
      }
    }
  }

    /*
  for (UInt n=0; n < this->active_impactor_nodes->getSize(); ++n) {
    UInt current_node = this->active_impactor_nodes->at(n);
    
    for (UInt i=0; i < spatial_dimension; ++i) {
      Int direction = this->master_normals->at(n,i);
      Real force = residual_val[current_node * spatial_dimension + i];
      if(force * direction > 0.) {
	bound_val[current_node * spatial_dimension + i] = false;
	this->active_impactor_nodes->erase(n);
	this->master_normals->erase(n);
	this->node_is_sticking->erase(n);
	n--;
	break;
      }
    }
    }*/
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addFriction() {
  AKANTU_DEBUG_IN();
  
  Real friction_coef = 0.3; // temp solution until friction coefficient better defined 
  
  const Real tolerance = std::numeric_limits<Real>::epsilon();
  
  Real * residual_val = this->model.getResidual().values;
  Real * velocity_val = this->model.getVelocity().values;
  
  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    // find index of master surface in impactors_information 
    Int master_index = -1;
    for (UInt i=0; i < this->impactors_information.size(); ++i) {
      if (this->impactors_information.at(i)->master_id == master) {
	master_index = i;
	break;
      }
    }
    AKANTU_DEBUG_ASSERT(master_index != -1, "No impactor information object for master" << master << "in impactors_information vector that is ");

    ImpactorInformationPerMaster * impactor_info = this->impactors_information.at(master_index);
    
    UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
    Int * direction_val = impactor_info->master_normals->values;
    bool * node_is_sticking_val = impactor_info->node_is_sticking->values;
    
    for (UInt n=0; n < impactor_info->active_impactor_nodes->getSize(); ++n) {
      UInt current_node = active_impactor_nodes_val[n];
      Real normal_contact_force = 0.;
      Real friction_force = 0.;
      
      // find friction force mu * normal force
      for (UInt i=0; i < spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] != 0) {
	  normal_contact_force = fabs(residual_val[current_node * this->spatial_dimension + i]);
	  friction_force = friction_coef * normal_contact_force;
	}
      }
      
      // find length of the residual projected to the frictional plane
      Real projected_residual = 0.;
      Real projected_velocity_magnitude = 0.;
      for (UInt i=0; i < this->spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] == 0) {
	  projected_residual += residual_val[current_node * this->spatial_dimension + i] * 
	                        residual_val[current_node * this->spatial_dimension + i];
	  projected_velocity_magnitude += velocity_val[current_node * this->spatial_dimension + i] *
	                                  velocity_val[current_node * this->spatial_dimension + i];
	}
      }
      projected_residual = sqrt(projected_residual);
      projected_velocity_magnitude = sqrt(projected_velocity_magnitude);
      
      // if it is a sticking node, check if it starts moving
      if(node_is_sticking_val[n*2+1]) {
	// node starts sliding
	if(projected_residual > friction_force)
	  node_is_sticking_val[n*2+1] = false;
	// node continues to stick and its friction force is equal the resiual
	else
	  friction_force = projected_residual;
      }
      
      // compute vector of length one in direction of projected residual
      Real * given_direction = NULL;
      Real given_length = 0.;
      if(node_is_sticking_val[n*2]) {
	given_direction = &residual_val[0];
	given_length = projected_residual;
      }
      else {
	given_direction = &velocity_val[0];
	given_length = projected_velocity_magnitude;
      }
      // if no tangential direction -> no friction force
      if(given_length < tolerance)
	continue;
      
      Real friction_direction[3];
      for (UInt i=0; i < this->spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] == 0)
	  friction_direction[i] = given_direction[current_node * this->spatial_dimension + i] / given_length;
	else
	  friction_direction[i] = 0.;
      }
      
      // add friction force to residual
      for (UInt i=0; i < this->spatial_dimension; ++i) 
	residual_val[current_node * this->spatial_dimension + i] -= friction_force * friction_direction[i];
      
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addSticking() {
  AKANTU_DEBUG_IN();
  
  Real * velocity_val = this->model.getVelocity().values;
  Real * acceleration_val = this->model.getAcceleration().values;
  const Real time_step = this->model.getTimeStep();

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    // find index of master surface in impactors_information 
    Int master_index = -1;
    for (UInt i=0; i < this->impactors_information.size(); ++i) {
      if (this->impactors_information.at(i)->master_id == master) {
	master_index = i;
	break;
      }
    }
    AKANTU_DEBUG_ASSERT(master_index != -1, "No impactor information object for master" << master << "in impactors_information vector that is ");

    ImpactorInformationPerMaster * impactor_info = this->impactors_information.at(master_index);
    
    UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
    Int * direction_val = impactor_info->master_normals->values;
    bool * node_is_sticking_val = impactor_info->node_is_sticking->values;

    for (UInt n=0; n < impactor_info->active_impactor_nodes->getSize(); ++n) {
      UInt current_node = active_impactor_nodes_val[n];
	
      if(!node_is_sticking_val[n*2]) {
	// compute scalar product of projected velocities
	Real scalar_prod_velocity = 0.;
	for (UInt i=0; i < this->spatial_dimension; ++i) {
	  if(direction_val[n * this->spatial_dimension + i] == 0) {
	    Real current_velocity = velocity_val[current_node * this->spatial_dimension + i];
	    Real estimated_velocity = current_velocity + time_step * acceleration_val[current_node * this->spatial_dimension + i];
	    scalar_prod_velocity += current_velocity * estimated_velocity;
	  }
	}
	// if scalar product <= 0, it has to be stick
	if(scalar_prod_velocity <= 0) {
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    if(direction_val[n * this->spatial_dimension + i] == 0) {
	      velocity_val[current_node * this->spatial_dimension + i]     = 0.;
	      acceleration_val[current_node * this->spatial_dimension + i] = 0.;
	    }
	  }
	  node_is_sticking_val[n*2]   = true;
	  node_is_sticking_val[n*2+1] = true;
	}
      }
    
      // for node that left sticking state set all sicking variables to false
      if(!node_is_sticking_val[n*2+1])
	node_is_sticking_val[n*2] = false;
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
