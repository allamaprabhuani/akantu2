/**
 * @file   contact_search.cc
 * @author David Kammer <kammer@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Fri Oct  8 11:46:34 2010
 *
 * @brief
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
#include "contact_search.hh"
#include "contact.hh"
#include "contact_neighbor_structure.hh"
#include "regular_grid_neighbor_structure.hh"
#include "grid_2d_neighbor_structure.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
PenetrationList::PenetrationList() : penetrating_nodes(Vector<UInt>(0, 1, "penetrating_nodes")) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PenetrationList::~PenetrationList() {
  AKANTU_DEBUG_IN();

  // for (UInt i = 0; i < _max_element_type; ++i) {
  //   if(penetrated_facets_offset[i]) delete penetrated_facets_offset[i];
  //   if(penetrated_facets       [i]) delete penetrated_facets       [i];
  //   if(facets_normals          [i]) delete facets_normals          [i];
  //   if(gaps                    [i]) delete gaps                    [i];
  //   if(projected_positions     [i]) delete projected_positions     [i];
  // }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactSearch::ContactSearch(Contact & contact,
			     const ContactNeighborStructureType & neighbors_structure_type,
			     const ContactSearchType & type,
			     const ContactSearchID & id) :
  id(id), contact(contact), neighbors_structure_type(neighbors_structure_type),
  type(type) {
  AKANTU_DEBUG_IN();

  contact.setContactSearch(*this);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ContactSearch::~ContactSearch() {
  AKANTU_DEBUG_IN();

  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    delete it->second;
  }

  neighbors_structure.clear();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::initSearch() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::initNeighborStructure() {
  AKANTU_DEBUG_IN();

  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    it->second->initNeighborStructure();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::initNeighborStructure(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  std::map<Surface, ContactNeighborStructure *>::iterator it;
  for (it = neighbors_structure.begin(); it != neighbors_structure.end(); ++it) {
    if(it->first == master_surface)
      it->second->initNeighborStructure();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::addMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) == neighbors_structure.end(),
		      "Master surface already registered in the search object " << id);

  ContactNeighborStructure * tmp_neighbors_structure = NULL;

  std::stringstream sstr;
  sstr << id << ":contact_neighbor_structure:" << neighbors_structure_type << ":" << master_surface;

  switch(neighbors_structure_type) {
  case _cnst_regular_grid : {
    Mesh & mesh = contact.getModel().getFEM().getMesh();
    if (mesh.getSpatialDimension() == 2) {
      tmp_neighbors_structure = new RegularGridNeighborStructure<2>(*this, master_surface, neighbors_structure_type, sstr.str());
    }
    else if(mesh.getSpatialDimension() == 3) {
      tmp_neighbors_structure = new RegularGridNeighborStructure<3>(*this, master_surface, neighbors_structure_type, sstr.str());
    }
    else
      AKANTU_DEBUG_ERROR("RegularGridNeighborStructure does not exist for dimension: "
			 << mesh.getSpatialDimension());
    break;
  }
  case _cnst_2d_grid : {
    tmp_neighbors_structure = new Grid2dNeighborStructure(*this, master_surface, neighbors_structure_type, sstr.str());
    break;
  }
  case _cnst_not_defined :
    //    tmp_neighbors_structure = new ContactNeighborStructureGrid2d(this, master_surface, sstr.str());
    AKANTU_DEBUG_ERROR("Not a valid neighbors structure type : " << neighbors_structure_type);
    break;
  }

  neighbors_structure[master_surface] = tmp_neighbors_structure;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ContactSearch::removeMasterSurface(const Surface & master_surface) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  delete neighbors_structure[master_surface];
  neighbors_structure.erase(neighbors_structure.find(master_surface));

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void ContactSearch::updateStructure(const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  neighbors_structure[master_surface]->update();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
bool ContactSearch::checkIfUpdateStructureNeeded(const Surface & master_surface) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(neighbors_structure.find(master_surface) != neighbors_structure.end(),
		      "Master surface not registered in the search object " << id);

  bool check = neighbors_structure[master_surface]->check();

  AKANTU_DEBUG_OUT();
  return check;
}

/* -------------------------------------------------------------------------- */
const ContactNeighborStructure & ContactSearch::getContactNeighborStructure(const Surface & master_surface) const {
  AKANTU_DEBUG_IN();
  std::map<Surface, ContactNeighborStructure *>::const_iterator it = neighbors_structure.find(master_surface);
  AKANTU_DEBUG_ASSERT(it != neighbors_structure.end(), "Master surface not registred in contact search.");

  AKANTU_DEBUG_OUT();
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
void ContactSearch::computeMaxIncrement(Real * max_increment) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = contact.getModel().getFEM().getMesh().getSpatialDimension();
  UInt nb_surfaces = contact.getModel().getFEM().getMesh().getNbSurfaces();
  Real * current_increment = contact.getModel().getIncrement().values;

  /// initialize max table with zeros
  for(UInt dim = 0; dim < spatial_dimension; ++dim)
    max_increment[dim] = 0.0;

  // get the nodes that are on the surfaces
  UInt * surface_to_nodes_offset = contact.getSurfaceToNodesOffset().values;
  UInt * surface_to_nodes        = contact.getSurfaceToNodes().values;

  /// find maximal increment of surface nodes in all directions
  for(UInt surf = 0; surf < nb_surfaces; ++surf) {
    UInt min_surf_offset = surface_to_nodes_offset[surf];
    UInt max_surf_offset = surface_to_nodes_offset[surf + 1];
    for(UInt n = min_surf_offset; n < max_surf_offset; ++n) {
      UInt cur_node = surface_to_nodes[n];
      for(UInt dim = 0; dim < spatial_dimension; ++dim) {
  	Real cur_increment = current_increment[cur_node * spatial_dimension + dim];
  	max_increment[dim] = std::max(max_increment[dim], fabs(cur_increment));
      }
    }
  }

  AKANTU_DEBUG_OUT();
}




__END_AKANTU__
