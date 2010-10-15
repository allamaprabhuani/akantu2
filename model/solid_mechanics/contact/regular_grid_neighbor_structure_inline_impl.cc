/**
 * @file   regular_grid_neighbor_structure_inline_impl.cc
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Oct 11 17:39:51 2010
 *
 * @brief Implementation of inline functions of the regular grid
 * neighbor structure class
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */
inline void setGridSpacing(Real spacing, UInt component) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");

  grid_spacing[component] = spacing;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline Real getGridSpacing(UInt component) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");
  
  AKANTU_DEBUG_OUT();
  return grid_spacing[component];
}

