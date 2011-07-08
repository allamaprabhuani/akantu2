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

  this->impactor_surfaces = new std::vector<Surface>(0);

  this->active_impactor_nodes = new Vector<UInt>(0,1);
  //this->master_element_offset = new Vector<UInt>(0,1);
  this->master_element_type = new std::vector<ElementType>(0);
  this->master_normals = new Vector<Real>(0, spatial_dimension);
  this->node_is_sticking = new Vector<bool>(0,2);
  this->friction_forces = new Vector<Real>(0,spatial_dimension);
  this->stick_positions = new Vector<Real>(0,spatial_dimension);
  this->residual_forces = new Vector<Real>(0,spatial_dimension);
  this->previous_velocities = new Vector<Real>(0,spatial_dimension);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactRigid::ImpactorInformationPerMaster::~ImpactorInformationPerMaster() {
  AKANTU_DEBUG_IN();

  delete this->impactor_surfaces;

  delete this->active_impactor_nodes;
  //  delete this->master_element_offset;
  delete this->master_element_type;
  delete this->master_normals;
  delete this->node_is_sticking;
  delete this->friction_forces;
  delete this->stick_positions;
  delete this->residual_forces;
  delete this->previous_velocities;

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

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  for (it_imp = this->impactors_information.begin();
       it_imp != this->impactors_information.end();
       ++it_imp) {
    delete it_imp->second;
  }

  this->impactors_information.clear();
  this->friction_coefficient.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  Contact::addMasterSurface(master_surface);

  ImpactorInformationPerMaster * impactor_info = new ImpactorInformationPerMaster(master_surface, this->spatial_dimension);
  this->impactors_information[master_surface] = impactor_info;

  //this->friction_coefficient[master_surface] = NULL;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addImpactorSurfaceToMasterSurface(const Surface & impactor_surface, const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master_surface);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master_surface << "couldn't be found and erased in impactors_information map");
  it_imp->second->impactor_surfaces->push_back(impactor_surface);

  /*
  for (UInt m=0; m < this->impactors_information.size(); ++m) {
    if (this->impactors_information.at(m)->master_id == master_surface) {
      this->impactors_information.at(m)->impactor_surfaces->push_back(impactor_surface);
      break;
    }
  }
  */

  ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
  it_fc = this->friction_coefficient.find(master_surface);
  if (it_fc != this->friction_coefficient.end())
    it_fc->second->addImpactorSurface(impactor_surface);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::removeMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  Contact::removeMasterSurface(master_surface);

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master_surface);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master_surface << "couldn't be found and erased in impactors_information map");
  delete it_imp->second;
  this->impactors_information.erase(it_imp);

  ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
  it_fc = this->friction_coefficient.find(master_surface);
  AKANTU_DEBUG_ASSERT(it_fc != this->friction_coefficient.end(),
		      "The master surface: " << master_surface << "couldn't be found and erased in friction_coefficient map");
  this->friction_coefficient.erase(it_fc);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::removeImpactorSurfaceFromMasterSurface(const Surface & impactor_surface, const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  // find in map the impactor information for the given master surface
  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master_surface);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master_surface << "couldn't be found and erased in impactors_information map");
  std::vector<Surface> * imp_surfaces = it_imp->second->impactor_surfaces;

  // find and delete the impactor surface
  std::vector<Surface>::iterator it_surf;
  it_surf = find(imp_surfaces->begin(), imp_surfaces->end(), impactor_surface);
  AKANTU_DEBUG_ASSERT(it_surf != imp_surfaces->end(), "Couldn't find impactor surface " << impactor_surface << " for the master surface " << master_surface << " and couldn't erase it");
  imp_surfaces->erase(it_surf);

  /*
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
  */

  ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
  it_fc = this->friction_coefficient.find(master_surface);
  if (it_fc != this->friction_coefficient.end())
    it_fc->second->removeImpactorSurface(impactor_surface);


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::solveContact() {
  AKANTU_DEBUG_IN();

  for(UInt m=0; m < master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    PenetrationList * penet_list = new PenetrationList();
    contact_search->findPenetration(master, *penet_list);

    solvePenetrationClosestProjection(master, *penet_list);
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
void ContactRigid::solvePenetrationClosestProjection(const Surface master,
						     const PenetrationList & penet_list) {
  AKANTU_DEBUG_IN();

  const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;

  /// find existing surface element types
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
      Real * gaps = penet_list.gaps(type, _not_ghost)->values;
      UInt offset_min = penet_list.penetrated_facets_offset(type, _not_ghost)->get(n);
      UInt offset_max = penet_list.penetrated_facets_offset(type, _not_ghost)->get(n+1);
      for (UInt f = offset_min; f < offset_max; ++f) {
	if(gaps[f] < min_gap) {
	  min_gap = gaps[f];
	  penetrated_type = type;
	  penetrated_facet_offset = f;
	}
      }
    }

    bool is_active = isAlreadyActiveImpactor(master, penet_list, n, penetrated_type, penetrated_facet_offset);
    if (!is_active) {
      // correct the position of the impactor
      projectImpactor(penet_list, n, penetrated_type, penetrated_facet_offset);
      lockImpactorNode(master, penet_list, n, penetrated_type, penetrated_facet_offset);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool ContactRigid::isAlreadyActiveImpactor(const Surface master,
					   const PenetrationList & penet_list,
					   const UInt impactor_index,
					   const ElementType facet_type,
					   const UInt facet_offset) {
  AKANTU_DEBUG_IN();

  bool is_active = false;

  UInt * penetrating_nodes = penet_list.penetrating_nodes.values;
  UInt impactor_node = penetrating_nodes[impactor_index];

  // find facet normal
  Real * facets_normals = penet_list.facets_normals(facet_type, _not_ghost)->values;
  Real * facet_normal = &facets_normals[facet_offset*spatial_dimension];
  Int normal[this->spatial_dimension];
  for(UInt i = 0; i < this->spatial_dimension; ++i)
    normal[i] = static_cast<Int>(floor(facet_normal[i] + 0.5));

  // check if this is already in the active impactor node list
  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master << "couldn't be found in impactors_information map");
  ImpactorInformationPerMaster * impactor_info = it_imp->second;
  Vector<UInt> * active_nodes = impactor_info->active_impactor_nodes;
  Vector<Real> * master_normals = impactor_info->master_normals;

  for (UInt i=0; i<active_nodes->getSize(); ++i) {
    if(active_nodes->at(i) == impactor_node) {
      UInt count = 0;
      for (UInt d=0; d<this->spatial_dimension; ++d) {
	if(master_normals->at(i,d) == normal[d])
	  count++;
      }
      if(count == this->spatial_dimension)
	is_active = true;
    }
  }

  AKANTU_DEBUG_OUT();
  return is_active;
}

/* -------------------------------------------------------------------------- */
void ContactRigid::projectImpactor(const PenetrationList & penet_list, const UInt impactor_index, const ElementType facet_type, const UInt facet_offset) {

  AKANTU_DEBUG_IN();

  const bool increment_flag = model.getIncrementFlag();

  UInt * penetrating_nodes   = penet_list.penetrating_nodes.values;
  Real * facets_normals      = penet_list.facets_normals(facet_type, _not_ghost)->values;
  Real * gaps                = penet_list.gaps(facet_type, _not_ghost)->values;
  Real * projected_positions = penet_list.projected_positions(facet_type, _not_ghost)->values;

  Real * current_position = model.getCurrentPosition().values;
  Real * displacement     = model.getDisplacement().values;
  Real * increment = NULL;
  if(increment_flag)
    increment = model.getIncrement().values;

  UInt impactor_node = penetrating_nodes[impactor_index];

  for(UInt i=0; i < this->spatial_dimension; ++i) {
    current_position[impactor_node*spatial_dimension + i] = projected_positions[facet_offset*spatial_dimension + i];
    Real displacement_correction = gaps[facet_offset] * facets_normals[facet_offset*spatial_dimension + i];
    displacement[impactor_node*spatial_dimension + i] = displacement[impactor_node*spatial_dimension + i] - displacement_correction;
    if(increment_flag)
      increment   [impactor_node*spatial_dimension + i] = increment   [impactor_node*spatial_dimension + i] - displacement_correction;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::lockImpactorNode(const Surface master, const PenetrationList & penet_list, const UInt impactor_index, const ElementType facet_type, const UInt facet_offset) {
  AKANTU_DEBUG_IN();

  bool init_state_stick = true; // node is in sticking when it gets in contact

  UInt * penetrating_nodes = penet_list.penetrating_nodes.values;
  UInt impactor_node = penetrating_nodes[impactor_index];

  Real * facets_normals = penet_list.facets_normals(facet_type, _not_ghost)->values;
  Real * facet_normal = &facets_normals[facet_offset*spatial_dimension];
  Real normal[this->spatial_dimension];
  //Int * normal_val = &normal[0];

  bool * bound_val    = this->model.getBoundary().values;
  Real * position_val = this->model.getCurrentPosition().values;
  Real * veloc_val    = this->model.getVelocity().values;
  Real * accel_val    = this->model.getAcceleration().values;
  Real * residual_val = this->model.getResidual().values;

  for(UInt i = 0; i < this->spatial_dimension; ++i)
    normal[i] = floor(facet_normal[i] + 0.5);

  for(UInt i = 0; i < this->spatial_dimension; ++i) {
    UInt index = impactor_node * spatial_dimension + i;
    if(normal[i] != 0) {
      bound_val[index] = true;
      veloc_val[index] = 0.;
      accel_val[index] = 0.;
    }
    // if node is in initial stick its tangential velocity has to be zero
    if(init_state_stick) {
      veloc_val[index] = 0.;
      accel_val[index] = 0.;
    }
  }

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master << "couldn't be found in impactors_information map");
  ImpactorInformationPerMaster * impactor_info = it_imp->second;
  impactor_info->active_impactor_nodes->push_back(impactor_node);
  //  impactor_info->master_element_offset->push_back(facet_offset);
  impactor_info->master_element_type->push_back(facet_type);
  impactor_info->master_normals->push_back(normal);
  // initial value for stick state when getting in contact
  bool init_sticking[2];
  if (init_state_stick) {
    init_sticking[0] = true;
    init_sticking[1] = true;
  }
  else {
    init_sticking[0] = false;
    init_sticking[1] = false;
  }
  impactor_info->node_is_sticking->push_back(init_sticking);
  Real init_friction_force[this->spatial_dimension];
  for(UInt i=0; i<this->spatial_dimension; ++i) {
    init_friction_force[i] = 0.;
  }
  impactor_info->friction_forces->push_back(init_friction_force);
  impactor_info->stick_positions->push_back(&(position_val[impactor_node*this->spatial_dimension]));
  impactor_info->residual_forces->push_back(&(residual_val[impactor_node*this->spatial_dimension]));
  impactor_info->previous_velocities->push_back(&(veloc_val[impactor_node*this->spatial_dimension]));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::avoidAdhesion() {
  AKANTU_DEBUG_IN();

  Real * residual_val = this->model.getResidual().values;
  bool * bound_val    = this->model.getBoundary().values;

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
    it_imp = this->impactors_information.find(master);
    AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
			"The master surface: " << master << "couldn't be found in impactors_information map");
    ImpactorInformationPerMaster * impactor_info = it_imp->second;

    for (UInt n=0; n < impactor_info->active_impactor_nodes->getSize(); ++n) {
      UInt current_node = impactor_info->active_impactor_nodes->at(n);

      for (UInt i=0; i < spatial_dimension; ++i) {
	Int direction = impactor_info->master_normals->at(n,i);
	Real force = residual_val[current_node * spatial_dimension + i];
	if(force * direction > 0.) {
	  bound_val[current_node * spatial_dimension + i] = false;
	  impactor_info->active_impactor_nodes->erase(n);
	  //	  impactor_info->master_element_offset->erase(n);
	  impactor_info->master_element_type->erase(impactor_info->master_element_type->begin()+n);
	  impactor_info->master_normals->erase(n);
	  impactor_info->node_is_sticking->erase(n);
	  impactor_info->friction_forces->erase(n);
	  impactor_info->stick_positions->erase(n);
	  impactor_info->residual_forces->erase(n);
	  impactor_info->previous_velocities->erase(n);
	  n--;
	  break;
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::frictionPredictor() {
  AKANTU_DEBUG_IN();

  const Real null_tolerance = std::numeric_limits<Real>::epsilon() * 2.;

  Real * residual_val     = this->model.getResidual().values;
  Real * acceleration_val = this->model.getAcceleration().values;
  Real * velocity_val     = this->model.getVelocity().values;

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
    it_imp = this->impactors_information.find(master);
    AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
			"The master surface: " << master << "couldn't be found in impactors_information map");
    ImpactorInformationPerMaster * impactor_info = it_imp->second;

    ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
    it_fc = this->friction_coefficient.find(master);
    AKANTU_DEBUG_ASSERT(it_fc != this->friction_coefficient.end(),
			"The master surface: " << master << "couldn't be found in friction_coefficient map");
    FrictionCoefficient * fric_coef = it_fc->second;
    AKANTU_DEBUG_ASSERT(fric_coef != NULL, "There is no friction coefficient defined for master surface " << master);

    UInt nb_active_impactor_nodes = impactor_info->active_impactor_nodes->getSize();

    // compute the friction coefficient for each active impactor node
    Vector<Real> friction_coefficient_values(nb_active_impactor_nodes, 1);
    Real * friction_coefficient_values_p = friction_coefficient_values.values;
    fric_coef->computeFrictionCoefficient(friction_coefficient_values);

    UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
    Real * direction_val = impactor_info->master_normals->values;
    bool * node_is_sticking_val = impactor_info->node_is_sticking->values;
    Real * friction_forces_val = impactor_info->friction_forces->values;
    Real * residual_forces_val = impactor_info->residual_forces->values;
    Real * previous_velocities_val = impactor_info->previous_velocities->values;

    for (UInt n=0; n < nb_active_impactor_nodes; ++n) {
      UInt current_node = active_impactor_nodes_val[n];

      // save original residual
      for(UInt i=0; i < spatial_dimension; ++i) {
	residual_forces_val[n*spatial_dimension+i] = residual_val[current_node*spatial_dimension+i];
      }

      // find friction force mu * normal force
      Real normal_contact_force = 0.;
      Real friction_force = 0.;
      for (UInt i=0; i < spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] != 0) {
	  normal_contact_force = fabs(residual_val[current_node * this->spatial_dimension + i]);
	  friction_force = friction_coefficient_values_p[n] * normal_contact_force;
	}
      }

      // if node is sticking, friction force acts against the residual
      if (node_is_sticking_val[n*2]) {
	// compute length of projected residual
	Real projected_residual = 0.;
	for (UInt i=0; i < this->spatial_dimension; ++i) {
	  if(direction_val[n * this->spatial_dimension + i] == 0) {
	    projected_residual += residual_val[current_node * this->spatial_dimension + i] *
				  residual_val[current_node * this->spatial_dimension + i];
	  }
	}
	projected_residual = sqrt(projected_residual);

	// friction is weaker than residual -> start sliding
	if (friction_force < projected_residual) {
	  node_is_sticking_val[n*2]   = false;
	  node_is_sticking_val[n*2+1] = false;

	  // compute friction vector
	  Real friction_direction[3];
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    if(direction_val[n * this->spatial_dimension + i] == 0) {
	      friction_direction[i] = residual_val[current_node * this->spatial_dimension + i];
	      friction_direction[i] /= projected_residual;
	    }
	    else
	      friction_direction[i] = 0.;
	    friction_direction[i] *= -1.;
	  }

	  // add friction force to residual
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    Real directional_fric_force = friction_force * friction_direction[i];
	    friction_forces_val[n*this->spatial_dimension + i] = directional_fric_force;
	    residual_val[current_node * this->spatial_dimension + i] += directional_fric_force;
	  }
	}

	// friction is stronger than residual -> stay sticking
	else {
	  // put residual in friction force and set residual to zero
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    if(direction_val[n * this->spatial_dimension + i] == 0) {
	      friction_forces_val[n*spatial_dimension + i] = residual_val[current_node * this->spatial_dimension + i];
	      residual_val[current_node * this->spatial_dimension + i] = 0.;
	    }
	    else
	      friction_forces_val[n*spatial_dimension + i]  = 0.;
	    friction_forces_val[n*spatial_dimension + i] *= -1.;
	  }
	}
      } // end if stick


      // if node is sliding friction force acts against velocity and
      // in case there is no velocity it acts against the residual
      else {
	Real projected_residual = 0.;
	Real projected_velocity_magnitude = 0.;
	for (UInt i=0; i < this->spatial_dimension; ++i) {
	  if(direction_val[n * this->spatial_dimension + i] == 0) {
	    projected_residual += residual_val[current_node * this->spatial_dimension + i] *
				  residual_val[current_node * this->spatial_dimension + i];
	    projected_velocity_magnitude += previous_velocities_val[n * this->spatial_dimension + i] *
					    previous_velocities_val[n * this->spatial_dimension + i];
	  }
	}
	projected_residual = sqrt(projected_residual);
	projected_velocity_magnitude = sqrt(projected_velocity_magnitude);

	// if no projected velocity nor friction force stronger than projected residual, this is stick
	if ((projected_velocity_magnitude < null_tolerance) && !(projected_residual > friction_force)) {
	  // put residual in friction force and set residual to zero
	  // set projected velocity and acceleration to zero
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    if(direction_val[n * this->spatial_dimension + i] == 0) {
	      friction_forces_val[n*spatial_dimension + i] = residual_val[current_node * this->spatial_dimension + i];
	      residual_val    [current_node * this->spatial_dimension + i] = 0.;
	      velocity_val    [current_node * this->spatial_dimension + i] = 0.;
	      acceleration_val[current_node * this->spatial_dimension + i] = 0.;
	    }
	    else
	      friction_forces_val[n*spatial_dimension + i] = 0.;
	    friction_forces_val[n*spatial_dimension + i] *= -1.;
	  }
	  node_is_sticking_val[n*2]   = true;
	  node_is_sticking_val[n*2+1] = true;
	}

	// node is really sliding
	else {
	  UInt index = 0;
	  Real * direction_vector;
	  Real direction_length;
	  // if only velocity is zero (because slider just turned direction)
	  if (projected_velocity_magnitude < null_tolerance) {
	    index = current_node;
	    direction_vector = residual_val;
	    direction_length = projected_residual;
	  }
	  // there is velocity
	  else {
	    index = n;
	    direction_vector = previous_velocities_val;
	    direction_length = projected_velocity_magnitude;
	  }

	  // compute the friction direction
	  Real friction_direction[3];
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    if(direction_val[n * this->spatial_dimension + i] == 0) {
	      friction_direction[i] = direction_vector[index * this->spatial_dimension + i];
	      friction_direction[i] /= direction_length;
	    }
	    else
	      friction_direction[i] = 0.;
	    friction_direction[i] *= -1.;
	  }

	  // add friction force to residual
	  for (UInt i=0; i < this->spatial_dimension; ++i) {
	    Real directional_fric_force = friction_force * friction_direction[i];
	    friction_forces_val[n*this->spatial_dimension + i] = directional_fric_force;
	    residual_val[current_node * this->spatial_dimension + i] += directional_fric_force;
	  }
	}
      } // end sliding
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::frictionCorrector() {
  AKANTU_DEBUG_IN();

  const Real tolerance = std::numeric_limits<Real>::epsilon() * 100.;

  Real * residual_val     = this->model.getResidual().values;
  Real * acceleration_val = this->model.getAcceleration().values;
  Real * velocity_val     = this->model.getVelocity().values;

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
    it_imp = this->impactors_information.find(master);
    AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
			"The master surface: " << master << "couldn't be found in impactors_information map");
    ImpactorInformationPerMaster * impactor_info = it_imp->second;

    UInt nb_active_impactor_nodes = impactor_info->active_impactor_nodes->getSize();

    UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
    Real * direction_val = impactor_info->master_normals->values;
    bool * node_is_sticking_val = impactor_info->node_is_sticking->values;
    Real * friction_forces_val = impactor_info->friction_forces->values;
    Real * residual_forces_val = impactor_info->residual_forces->values;
    Real * previous_velocities_val = impactor_info->previous_velocities->values;

    for (UInt n=0; n < nb_active_impactor_nodes; ++n) {
      UInt current_node = active_impactor_nodes_val[n];

      // check only sliding impactor nodes
      if (node_is_sticking_val[n*2])
	continue;

      // compute scalar product of previous velocity with current velocity
      Real dot_prod_velocities = 0.;
      Real current_proj_vel_mag = 0.;
      for (UInt i=0; i < this->spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] == 0) {
	  dot_prod_velocities += velocity_val[current_node * this->spatial_dimension + i] *
				 previous_velocities_val[n * this->spatial_dimension + i];
	  current_proj_vel_mag += velocity_val[current_node * this->spatial_dimension + i] *
				  velocity_val[current_node * this->spatial_dimension + i];
	  }
      }
      current_proj_vel_mag = sqrt(current_proj_vel_mag);

      // if velocity has changed sign or has become very small we get into stick
      if ((dot_prod_velocities < 0.) || (current_proj_vel_mag < tolerance)) {
	// compute scalar product of residual and current velocity
	Real dot_prod_vel_res = 0.;
	Real proj_friction_mag = 0.;
	Real proj_residual_mag = 0.;
	for (UInt i=0; i < this->spatial_dimension; ++i) {
	  if(direction_val[n * this->spatial_dimension + i] == 0) {
	    dot_prod_vel_res += residual_forces_val[n * this->spatial_dimension + i] *
				previous_velocities_val[n * this->spatial_dimension + i];
	    proj_friction_mag += friction_forces_val[n * this->spatial_dimension + i] *
				 friction_forces_val[n * this->spatial_dimension + i];
	    proj_residual_mag += residual_forces_val[n * this->spatial_dimension + i] *
				 residual_forces_val[n * this->spatial_dimension + i];
	  }
	}
	proj_friction_mag = sqrt(proj_friction_mag);
	proj_residual_mag = sqrt(proj_residual_mag);

	for (UInt i=0; i < this->spatial_dimension; ++i) {
	  if(direction_val[n * this->spatial_dimension + i] == 0) {
	    //friction_forces_val[n*spatial_dimension + i] = residual_forces_val[n*spatial_dimension + i];
	    residual_forces_val[n * this->spatial_dimension + i] =0.;
	    residual_val    [current_node * this->spatial_dimension + i] = 0.;
	    velocity_val    [current_node * this->spatial_dimension + i] = 0.;
	    acceleration_val[current_node * this->spatial_dimension + i] = 0.;
	  }
	  else
	    friction_forces_val[n*spatial_dimension + i] = 0.;
	  friction_forces_val[n*spatial_dimension + i] *= -1.;
	}
	node_is_sticking_val[n*2]   = true;
	node_is_sticking_val[n*2+1] = true;
      }
    }

    // save current velocities
    for (UInt n=0; n < nb_active_impactor_nodes; ++n) {
      UInt current_node = active_impactor_nodes_val[n];
      for(UInt i=0; i < spatial_dimension; ++i) {
	previous_velocities_val[n*spatial_dimension+i] = velocity_val[current_node*spatial_dimension+i];
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::addRegularizedFriction(const Real & regularizer) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(this->spatial_dimension == 2,
			"addRegularizedFriction is implemented only for 2D");

  Real * residual_val = this->model.getResidual().values;
  Real * position_val = this->model.getCurrentPosition().values;

  for(UInt m=0; m < this->master_surfaces.size(); ++m) {
    Surface master = this->master_surfaces.at(m);

    // find the impactors information for this master surface
    ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
    it_imp = this->impactors_information.find(master);
    AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
			"The master surface: " << master << "couldn't be found in impactors_information map");
    ImpactorInformationPerMaster * impactor_info = it_imp->second;

    // find the friction coefficient object for this master
    ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
    it_fc = this->friction_coefficient.find(master);
    AKANTU_DEBUG_ASSERT(it_fc != this->friction_coefficient.end(),
			"The master surface: " << master << "couldn't be found in friction_coefficient map");
    FrictionCoefficient * fric_coef = it_fc->second;
    AKANTU_DEBUG_ASSERT(fric_coef != NULL, "There is no friction coefficient defined for master surface " << master);

    UInt nb_active_impactor_nodes = impactor_info->active_impactor_nodes->getSize();

    // compute the friction coefficient for each active impactor node
    Vector<Real> friction_coefficient_values(nb_active_impactor_nodes, 1);
    Real * friction_coefficient_values_p = friction_coefficient_values.values;
    fric_coef->computeFrictionCoefficient(friction_coefficient_values);

    UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
    Real * direction_val = impactor_info->master_normals->values;
    bool * node_is_sticking_val = impactor_info->node_is_sticking->values;
    Real * friction_forces_val = impactor_info->friction_forces->values;
    Real * stick_positions_val = impactor_info->stick_positions->values;

    for (UInt n=0; n < nb_active_impactor_nodes; ++n) {
      UInt current_node = active_impactor_nodes_val[n];
      Real normal_contact_force = 0.;
      Real friction_force = 0.;

      // find friction force mu * normal force
      for (UInt i=0; i<spatial_dimension; ++i) {
	if(direction_val[n * this->spatial_dimension + i] != 0) {
	  normal_contact_force = fabs(residual_val[current_node * this->spatial_dimension + i]);
	  friction_force = friction_coefficient_values_p[n] * normal_contact_force;
	}
      }

      // compute F_(i+1)trial
      Real f_trial = Math::distance_2d(&(position_val[current_node * this->spatial_dimension]),
				       &(stick_positions_val[n * this->spatial_dimension]));
      f_trial *= regularizer;

      Real delta_s_vector[this->spatial_dimension];
      for (UInt i=0; i<spatial_dimension; ++i) {
	delta_s_vector[i] = position_val[current_node * this->spatial_dimension + i]
			  - stick_positions_val[n * this->spatial_dimension + i];
	delta_s_vector[i] *= -1.;
      }

      // if node is on its own stick position no need to compute friction force
      // this avoids nan on normalize2
      if(Math::norm2(delta_s_vector) < Math::tolerance) {
	node_is_sticking_val[n*2] = true; node_is_sticking_val[n*2+1] = true;
	continue;
      }
      Math::normalize2(delta_s_vector);

      // check if node is still sticking
      if (f_trial <= friction_force) {
	friction_force = f_trial;
	// sticking position stays unchanged
	node_is_sticking_val[n*2] = true; node_is_sticking_val[n*2+1] = true;
      }
      else {
	node_is_sticking_val[n*2] = false; node_is_sticking_val[n*2+1] = false;
	Real delta_s = std::abs(f_trial - friction_force) / regularizer;
	// friction force stays unchanged
	for (UInt i=0; i<spatial_dimension; ++i) {
	  stick_positions_val[n * this->spatial_dimension + i] -= delta_s * delta_s_vector[i];
	}
      }

      // add friction force to residual
      for (UInt i=0; i < this->spatial_dimension; ++i) {
	Real directional_fric_force = friction_force * delta_s_vector[i];
	friction_forces_val[n*this->spatial_dimension + i] = directional_fric_force;
	residual_val[current_node * this->spatial_dimension + i] += directional_fric_force;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::setStickPositionsToCurrentPositions(const Surface master) {
  AKANTU_DEBUG_IN();

  Real * position_val = this->model.getCurrentPosition().values;

  // find the impactors information for this master surface
  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master << "couldn't be found in impactors_information map");
  ImpactorInformationPerMaster * impactor_info = it_imp->second;

  UInt nb_active_impactor_nodes = impactor_info->active_impactor_nodes->getSize();
  UInt * active_impactor_nodes_val = impactor_info->active_impactor_nodes->values;
  Real * stick_positions_val = impactor_info->stick_positions->values;

  for (UInt n=0; n < nb_active_impactor_nodes; ++n) {
    UInt current_node = active_impactor_nodes_val[n];

    // find friction force mu * normal force
    for (UInt i=0; i<spatial_dimension; ++i) {
      stick_positions_val[n*this->spatial_dimension + i] = position_val[current_node*this->spatial_dimension + i];
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::setFrictionCoefficient(FrictionCoefficient * fric_coefficient) {
  AKANTU_DEBUG_IN();

  Surface master = fric_coefficient->getMasterSurface();

  ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
  it_fc = this->friction_coefficient.find(master);
  AKANTU_DEBUG_ASSERT(it_fc == this->friction_coefficient.end(),
		      "There is already a friction coefficient for master surface: " << master);
  this->friction_coefficient[master] = fric_coefficient;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::removeFrictionCoefficient(const Surface master) {
  AKANTU_DEBUG_IN();

  ContactRigid::SurfaceToFricCoefMap::iterator it_fc;
  it_fc = this->friction_coefficient.find(master);
  AKANTU_DEBUG_ASSERT(it_fc != this->friction_coefficient.end(),
		      "There is no friction coefficient for master surface: " << master);
  this->friction_coefficient.erase(it_fc);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::getRestartInformation(std::map<std::string, VectorBase* > & map_to_fill,
					 Surface master) {
  AKANTU_DEBUG_IN();

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master << "couldn't be found in impactors_information map");
  ImpactorInformationPerMaster * impactor_info = it_imp->second;
  UInt nb_active_nodes = impactor_info->active_impactor_nodes->getSize();

  // create vectors
  Vector<bool> * activity_of_nodes = new Vector<bool>(this->mesh.getNbNodes(), 1, false);
  Vector<ElementType> * element_type_of_nodes = new Vector<ElementType>(this->mesh.getNbNodes(), 1, _not_defined);
  Vector<Real> * normals_of_nodes = new Vector<Real>(this->mesh.getNbNodes(), this->spatial_dimension, 0.);
  Vector<bool> * sticking_of_nodes = new Vector<bool>(this->mesh.getNbNodes(), 2, false);
  Vector<Real> * friction_force_of_nodes = new Vector<Real>(this->mesh.getNbNodes(), this->spatial_dimension, 0.);
  Vector<Real> * stick_position_of_nodes = new Vector<Real>(this->mesh.getNbNodes(), this->spatial_dimension, 0.);
  Vector<Real> * residual_force_of_nodes = new Vector<Real>(this->mesh.getNbNodes(), this->spatial_dimension, 0.);
  Vector<Real> * previous_velocity_of_nodes = new Vector<Real>(this->mesh.getNbNodes(), this->spatial_dimension, 0.);

  UInt * active_nodes = impactor_info->active_impactor_nodes->values;
  ElementType * element_type = &(*impactor_info->master_element_type)[0];
  Real * master_normal = impactor_info->master_normals->values;
  bool * node_stick = impactor_info->node_is_sticking->values;
  Real * friction_force = impactor_info->friction_forces->values;
  Real * stick_position = impactor_info->stick_positions->values;
  Real * residual_force = impactor_info->residual_forces->values;
  Real * previous_velocity = impactor_info->previous_velocities->values;

  for (UInt i=0; i<nb_active_nodes; ++i, ++active_nodes, ++element_type) {
    (*activity_of_nodes)(*active_nodes) = true;
    (*element_type_of_nodes)(*active_nodes) = *element_type;
    for (UInt d=0; d<this->spatial_dimension; ++d, ++master_normal, ++friction_force, ++stick_position, ++residual_force, ++previous_velocity) {
      (*normals_of_nodes)(*active_nodes,d) = *master_normal;
      (*friction_force_of_nodes)(*active_nodes,d) = *friction_force;
      (*stick_position_of_nodes)(*active_nodes,d) = *stick_position;
      (*residual_force_of_nodes)(*active_nodes,d) = *residual_force;
      (*previous_velocity_of_nodes)(*active_nodes,d) = *previous_velocity;
    }
    for (UInt j=0; j<2; ++j, ++node_stick) {
      (*sticking_of_nodes)(*active_nodes,j) = *node_stick;
    }
  }

  map_to_fill["active_impactor_nodes"] = activity_of_nodes;
  map_to_fill["master_element_type"] = element_type_of_nodes;
  map_to_fill["master_normals"] = normals_of_nodes;
  map_to_fill["node_is_sticking"] = sticking_of_nodes;
  map_to_fill["friction_forces"] = friction_force_of_nodes;
  map_to_fill["stick_positions"] = stick_position_of_nodes;
  map_to_fill["residual_forces"] = residual_force_of_nodes;
  map_to_fill["previous_velocities"] = previous_velocity_of_nodes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigid::setRestartInformation(std::map<std::string, VectorBase* > & restart_map,
					 Surface master) {
  AKANTU_DEBUG_IN();

  ContactRigid::SurfaceToImpactInfoMap::iterator it_imp;
  it_imp = this->impactors_information.find(master);
  AKANTU_DEBUG_ASSERT(it_imp != this->impactors_information.end(),
		      "The master surface: " << master << "couldn't be found in impactors_information map");
  ImpactorInformationPerMaster * impactor_info = it_imp->second;

  std::map < std::string, VectorBase* >::iterator it;

  it = restart_map.find("active_impactor_nodes");
  Vector<bool> * ai_nodes = (Vector<bool> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry for active impactor nodes" << std::endl;
  }

  it = restart_map.find("master_element_type");
  Vector<ElementType> * et_nodes = (Vector<ElementType> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry master element type" << std::endl;
  }

  it = restart_map.find("master_normals");
  Vector<Real> * mn_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry for master normals" << std::endl;
  }

  it = restart_map.find("node_is_sticking");
  Vector<bool> * is_nodes = (Vector<bool> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry node is sticking" << std::endl;
  }

  it = restart_map.find("friction_forces");
  Vector<Real> * ff_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry friction forces" << std::endl;
  }

  it = restart_map.find("stick_positions");
  Vector<Real> * sp_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry stick positions" << std::endl;
  }

  it = restart_map.find("residual_forces");
  Vector<Real> * rf_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry residual forces" << std::endl;
  }

  it = restart_map.find("previous_velocities");
  Vector<Real> * pv_nodes = (Vector<Real> *)(it->second);
  if(it == restart_map.end()) {
    std::cout << "could not find map entry previous" << std::endl;
  }



  UInt nb_nodes = this->mesh.getNbNodes();
  for (UInt i=0; i<nb_nodes; ++i) {
    if ((*ai_nodes)(i)) {
      // get active_impactor_nodes and master_element_type information
      impactor_info->active_impactor_nodes->push_back(i);
      impactor_info->master_element_type->push_back((*et_nodes)(i));

      // get master_normals and friction_forces information
      Real normal[this->spatial_dimension];
      Real friction[this->spatial_dimension];
      Real position[this->spatial_dimension];
      Real residual[this->spatial_dimension];
      Real velocity[this->spatial_dimension];
      for (UInt d=0; d<this->spatial_dimension; ++d) {
	normal[d] = (*mn_nodes)(i,d);
	friction[d] = (*ff_nodes)(i,d);
	position[d] = (*sp_nodes)(i,d);
	residual[d] = (*rf_nodes)(i,d);
	velocity[d] = (*pv_nodes)(i,d);
      }
      impactor_info->master_normals->push_back(normal);
      impactor_info->friction_forces->push_back(friction);
      impactor_info->stick_positions->push_back(position);
      impactor_info->residual_forces->push_back(residual);
      impactor_info->previous_velocities->push_back(velocity);

      // get node_is_sticking information
      bool stick[2];
      for (UInt j=0; j<2; ++j) {
	stick[j] = (*is_nodes)(i,j);
      }
      impactor_info->node_is_sticking->push_back(stick);
    }
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__
