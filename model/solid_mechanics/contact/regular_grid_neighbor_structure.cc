/**
 * @file   regular_grid_neighbor_structure.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Oct 11 16:03:17 2010
 *
 * @brief  Specialization of the contact neighbor structure for regular grid
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
#include "regular_grid_neighbor_structure.hh"


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
RegularGridNeighborStructure::RegularGridNeighborStructure(const ContactSearch & contact_search,
							   const Surface & master_surface,
							   const ContactNeighborStructureID & id) :
  ContactNeighborStructure(contact_search, master_surface, id) {
  AKANTU_DEBUG_IN();
  spatial_dimension = contact_search.getContact().getModel().getFEM().getMesh().getSpatialDimension();
  grid_spacing[0] = 0.1;
  grid_spacing[1] = 0.1;
  grid_spacing[2] = 0.1;
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void RegularGridNeighborStructure::init() {
  AKANTU_DEBUG_IN();
  this->update();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void RegularGridNeighborStructure::update() {
  AKANTU_DEBUG_IN();
  
  

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */




__END_AKANTU__
