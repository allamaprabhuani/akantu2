/**
 * @file   ruina_slowness_fric_coef.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Mar  7 16:19:04 2011
 *
 * @brief  implementation of methods for slowness law of the theta state variable
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
#include "ruina_slowness_fric_coef.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
RuinaSlownessFricCoef<compute_analytic_solution>::RuinaSlownessFricCoef(ContactRigid & contact,
									const Surface & master_surface) : SimplifiedDieterichFricCoef(contact, master_surface), d_zero(1.) {
  AKANTU_DEBUG_IN();

  // find in map the impactor information for the given master surface  
  const ContactRigid::SurfaceToImpactInfoMap & imp_info = this->contact.getImpactorsInformation(); 
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it;
  it = imp_info.find(master_surface);
  AKANTU_DEBUG_ASSERT(it != imp_info.end(), 
		      "The master surface: " << master_surface << "couldn't be found in impactors_information map");
  std::vector<Surface> * imp_surfaces = it->second->impactor_surfaces;

  const Vector<UInt> surface_to_nodes_offset = this->contact.getSurfaceToNodesOffset();
  const Vector<UInt> surface_to_nodes = this->contact.getSurfaceToNodes();
  UInt * surface_to_nodes_offset_val = surface_to_nodes_offset.values;
  UInt * surface_to_nodes_val = surface_to_nodes.values;

  for(UInt s = 0; s < imp_surfaces->size(); ++s) {
    UInt surf = imp_surfaces->at(s);
    UInt min_surf_offset = surface_to_nodes_offset_val[surf];
    UInt max_surf_offset = surface_to_nodes_offset_val[surf+1];
    for (UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      this->node_to_index[surface_to_nodes_val[n]] = this->node_to_index.size();
    }
  }
  
  UInt node_to_index_size = this->node_to_index.size();
  this->are_active_impactor_nodes = new Vector<bool>(node_to_index_size, 1, false);
  this->were_active_impactor_nodes = new Vector<bool>(node_to_index_size, 1, false);
  this->previous_theta_state_variables = new Vector<Real>(node_to_index_size, 1, 0.);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
RuinaSlownessFricCoef<compute_analytic_solution>::~RuinaSlownessFricCoef() {
  AKANTU_DEBUG_IN();

  delete this->are_active_impactor_nodes;
  delete this->were_active_impactor_nodes;
  delete this->previous_theta_state_variables;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
void RuinaSlownessFricCoef<compute_analytic_solution>::initializeComputeFricCoef() {
  AKANTU_DEBUG_IN();
  
  SimplifiedDieterichFricCoef::initializeComputeFricCoef();

  // theta value for a node that was not in contact before
  Real init_theta = 0.;

  // find impactor_information for given master
  const ContactRigid::SurfaceToImpactInfoMap & imp_info = this->contact.getImpactorsInformation();
  ContactRigid::SurfaceToImpactInfoMap::const_iterator it;

  it = imp_info.find(this->master_surface);
  AKANTU_DEBUG_ASSERT(it != imp_info.end(), 
		      "Couldn't find impactor information object for master surface " << master_surface);
  ContactRigid::ImpactorInformationPerMaster * impactor_info = it->second;

  // find the current active impactor nodes
  Vector<UInt> * active_nodes = impactor_info->active_impactor_nodes;
  UInt * active_nodes_val = active_nodes->values;

  UInt nb_impactor_nodes = this->node_to_index.size();
  bool * are_active_val = this->are_active_impactor_nodes->values;
  bool * were_active_val = this->were_active_impactor_nodes->values;
  Real * previous_theta_val = this->previous_theta_state_variables->values;

  // loop over current active impactor nodes and compute the new theta value
  this->theta_state_variables->resize(active_nodes->getSize());
  Real * theta_state_variables_val = this->theta_state_variables->values;

  // get the sliding speed for all active impactor nodes
  Real * relative_sliding_velocities_val = this->relative_sliding_velocities->values;

  Real delta_t = this->contact.getModel().getTimeStep();

  std::map<UInt, UInt>::iterator index_it;
  
  for(UInt n = 0; n < active_nodes->getSize(); ++n) {
    UInt impactor = active_nodes_val[n];
    index_it = this->node_to_index.find(impactor);
    AKANTU_DEBUG_ASSERT(index_it != this->node_to_index.end(), "Could not find impactor node " << impactor << " in node_to_index map of friction coefficient object");
    UInt index = index_it->second;
    
    // if nodes was not in contact before set thetas to initial value
    if(!were_active_val[index]) {
      theta_state_variables_val[n] = init_theta;
      previous_theta_val[index] = init_theta;
    }
    else {
      Real theta;
      if (compute_analytic_solution)
	theta = computeAnalyticTheta(previous_theta_val[index], 
				     relative_sliding_velocities_val[n],
				     delta_t);
      else
	theta = computeImplicitTheta(previous_theta_val[index], 
				     relative_sliding_velocities_val[n],
				     delta_t);
      
      theta_state_variables_val[n] = theta;
      previous_theta_val[index] = theta;
    }
    are_active_val[index] = true;
  }

  // set to zero all thetas which were not computed at this time step
  Vector<bool> * tmp;
  tmp = this->were_active_impactor_nodes;
  this->were_active_impactor_nodes = this->are_active_impactor_nodes;
  this->are_active_impactor_nodes = tmp;
  memset(this->are_active_impactor_nodes->values, false, nb_impactor_nodes*sizeof(bool));
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
void RuinaSlownessFricCoef<compute_analytic_solution>::addImpactorSurface(const Surface & impactor_surface) {
  AKANTU_DEBUG_IN();

  UInt init_note_to_index_size = this->node_to_index.size();

  // find pointer to the end (after last element) of vectors
  bool * end_node_ptr = &(this->were_active_impactor_nodes->values[this->were_active_impactor_nodes->getSize()]);
  Real * end_theta_ptr = &(this->previous_theta_state_variables->values[this->previous_theta_state_variables->getSize()]);

  // add new impactor nodes to the map
  const Vector<UInt> surface_to_nodes_offset = this->contact.getSurfaceToNodesOffset();
  const Vector<UInt> surface_to_nodes = this->contact.getSurfaceToNodes();
  UInt * surface_to_nodes_offset_val = surface_to_nodes_offset.values;
  UInt * surface_to_nodes_val = surface_to_nodes.values;

  UInt min_surf_offset = surface_to_nodes_val[impactor_surface];
  UInt max_surf_offset = surface_to_nodes_val[impactor_surface+1];
  for (UInt n = min_surf_offset; n < max_surf_offset; ++n) {
    this->node_to_index[surface_to_nodes_val[n]] = this->node_to_index.size();
  }

  // resize and initialize the vectors
  UInt node_to_index_size = this->node_to_index.size();
  UInt size_diff = node_to_index_size - init_note_to_index_size;
  if (size_diff > 0) {
    this->were_active_impactor_nodes->resize(node_to_index_size);
    memset(end_node_ptr, false, size_diff*sizeof(bool));
    this->previous_theta_state_variables->resize(node_to_index_size);
    memset(end_theta_ptr, false, size_diff*sizeof(Real));
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
void RuinaSlownessFricCoef<compute_analytic_solution>::removeImpactorSurface(const Surface & impactor_surface) {
  AKANTU_DEBUG_IN();

  /* do not take out the impactor nodes, because we do not know if they belong only to the given impactor_surface or if they belong also to another impactor_surface */

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<bool compute_analytic_solution>
void RuinaSlownessFricCoef<compute_analytic_solution>::setParam(const std::string & key, const std::string & value) {
  AKANTU_DEBUG_IN();
  
  std::stringstream sstr(value);
  if(key == "d_zero") { sstr >> this->d_zero; }
  //else if(key == "compute_analytic_solution") { sstr >> this->compute_analytic_solution; }
  else { SimplifiedDieterichFricCoef::setParam(key, value); }
  
  AKANTU_DEBUG_OUT();
}

template class RuinaSlownessFricCoef<true>;
template class RuinaSlownessFricCoef<false>;

__END_AKANTU__
