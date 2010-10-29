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
template<UInt spatial_dimension>
inline void RegularGridNeighborStructure<spatial_dimension>::setGridSpacing(Real spacing, UInt component) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");

  grid_spacing[component] = spacing;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real RegularGridNeighborStructure<spatial_dimension>::getGridSpacing(UInt component) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");
  
  AKANTU_DEBUG_OUT();
  return grid_spacing[component];
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void RegularGridNeighborStructure<spatial_dimension>::setSecurityFactor(Real factor, UInt component) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");

  security_factor[component] = factor;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real RegularGridNeighborStructure<spatial_dimension>::getSecurityFactor(UInt component) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " << 
		      component << " is out of range (spatial dimension = " << 
		      spatial_dimension << ")");
  
  AKANTU_DEBUG_OUT();
  return security_factor[component];
}


