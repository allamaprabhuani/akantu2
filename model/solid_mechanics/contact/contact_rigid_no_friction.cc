/**
 * @file   contact_rigid_no_friction.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Tue Oct 26 17:15:08 2010
 *
 * @brief Specialization  of the  contact structure for  3d contact  in explicit
 * time scheme
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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
#include "contact_rigid_no_friction.hh"
#include "contact_search.hh"


__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
ContactRigidNoFriction::ContactRigidNoFriction(const SolidMechanicsModel & model,
					       const ContactType & type,
					       const ContactID & id,
					       const MemoryID & memory_id) :
  Contact(model, type, id, memory_id), spatial_dimension(model.getSpatialDimension()), mesh(model.getFEM().getMesh()) {
  AKANTU_DEBUG_IN();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactRigidNoFriction::~ContactRigidNoFriction() {
  AKANTU_DEBUG_IN();

  

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactRigidNoFriction::solveContact() {
  AKANTU_DEBUG_IN();

  for(UInt master=0; master < master_surfaces.size(); ++master) {
    PenetrationList * penet_list = new PenetrationList();
    contact_search->findPenetration(master_surfaces.at(master), *penet_list);
    solvePenetrationClosestProjection(*penet_list);
    delete penet_list;
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*void ContactRigidNoFriction::solvePenetration(const PenetrationList & penet_list) {
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
void ContactRigidNoFriction::solvePenetrationClosestProjection(const PenetrationList & penet_list) {
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
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactRigidNoFriction::projectImpactor(const PenetrationList & penet_list, const UInt impactor_index, const ElementType facet_type, const UInt facet_offset) {

  AKANTU_DEBUG_IN();
  
  const UInt dim = model.getSpatialDimension();
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

  for(UInt i=0; i < dim; ++i) {
    current_position[impactor_node*dim + i] = projected_positions[facet_offset*dim + i];
    Real displacement_correction = gaps[facet_offset] * facets_normals[facet_offset*dim + i];
    displacement[impactor_node*dim + i] = displacement[impactor_node*dim + i] - displacement_correction;
    if(increment_flag)
      increment   [impactor_node*dim + i] = increment   [impactor_node*dim + i] - displacement_correction;
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__
