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

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void RegularGridNeighborStructure<spatial_dimension>::setMaxIncrement(Real increment, UInt component) {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " <<
		      component << " is out of range (spatial dimension = " <<
		      spatial_dimension << ")");

  max_increment[component] = increment;
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline Real RegularGridNeighborStructure<spatial_dimension>::getMaxIncrement(UInt component) const {
  AKANTU_DEBUG_IN();
  AKANTU_DEBUG_ASSERT(component < spatial_dimension, "The component " <<
		      component << " is out of range (spatial dimension = " <<
		      spatial_dimension << ")");

  AKANTU_DEBUG_OUT();
  return max_increment[component];
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt RegularGridNeighborStructure<spatial_dimension>::computeCellNb(Int * directional_nb_cells,
									   Int * directional_cell) {

  AKANTU_DEBUG_IN();
  UInt cell_number = directional_cell[spatial_dimension - 1];
  for(Int dim = spatial_dimension - 2; dim >= 0; --dim) {
    cell_number *= directional_nb_cells[dim];
    cell_number += directional_cell[dim];
  }

  AKANTU_DEBUG_OUT();
  return cell_number;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt RegularGridNeighborStructure<spatial_dimension>::computeNeighborCells(UInt cell,
										  UInt * neighbors,
										  Int * directional_nb_cells) {
  AKANTU_DEBUG_IN();
  UInt nb_neighbors = 0;
  UInt max_spatial_dimension = 3; // to avoid warnings while compiling
  Int directional_cell[max_spatial_dimension];
  UInt global_cell_nb = cell;

  /// find the directional cell number
  //for(Int dir = spatial_dimension - 1; dir >= 0; --dir) {
  for(UInt dir = spatial_dimension - 1; dir != std::numeric_limits<UInt>::max(); --dir) {
    UInt factor = 1;
    for(UInt i = 0; i < dir; ++i) {
      factor *= directional_nb_cells[i];
    }
    directional_cell[dir] = std::floor(global_cell_nb / factor); // integer division !
    global_cell_nb -= directional_cell[dir] * factor;
  }

  /// compute neighbor cells
  Int neighbor_thickness = 1; // the number of neighbors for a given direction

  /// computation for 2D
  if(spatial_dimension == 2) {

    for(Int x = directional_cell[0] - neighbor_thickness; x <= directional_cell[0] + neighbor_thickness; ++x) {
      if(x < 0 || x >= directional_nb_cells[0]) continue; // border cell?

      for(Int y = directional_cell[1] - neighbor_thickness; y <= directional_cell[1] + neighbor_thickness; ++y) {
	if(y < 0 || y >= directional_nb_cells[1]) continue; // border cell?
	if(x == directional_cell[0] && y == directional_cell[1]) continue; // do only return neighbors not itself!

	Int neighbor_directional_cell[2];
	neighbor_directional_cell[0] = x;
	neighbor_directional_cell[1] = y;

	/// compute global cell index
	UInt neighbor_cell = computeCellNb(directional_nb_cells, neighbor_directional_cell);

	/// add the neighbor cell to the list
	neighbors[nb_neighbors++] = neighbor_cell;
      }
    }
  }
  /// computation for 3D
  else if(spatial_dimension == 3) {

    for(Int x = directional_cell[0] - neighbor_thickness; x <= directional_cell[0] + neighbor_thickness; ++x) {
      if(x < 0 || x >= directional_nb_cells[0]) continue; // border cell?

      for(Int y = directional_cell[1] - neighbor_thickness; y <= directional_cell[1] + neighbor_thickness; ++y) {
	if(y < 0 || y >= directional_nb_cells[1]) continue; // border cell?

	for(Int z = directional_cell[2] - neighbor_thickness; z <= directional_cell[2] + neighbor_thickness; ++z) {
	  if(z < 0 || z >= directional_nb_cells[2]) continue; // border cell?
	  if(x == directional_cell[0] && y == directional_cell[1] && z == directional_cell[2]) continue; // do only return neighbors not itself!

	  UInt local_spatial_dimension = 3; // to avoid warnings while compiling
	  Int neighbor_directional_cell[local_spatial_dimension];
	  neighbor_directional_cell[0] = x;
	  neighbor_directional_cell[1] = y;
	  neighbor_directional_cell[2] = z;

	  /// compute global cell index
	  UInt neighbor_cell = computeCellNb(directional_nb_cells, neighbor_directional_cell);

	  /// add the neighbor cell to the list
	  neighbors[nb_neighbors++] = neighbor_cell;
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
  return nb_neighbors;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void RegularGridNeighborStructure<spatial_dimension>::constructNeighborList() {

  AKANTU_DEBUG_IN();
  std::stringstream sstr; sstr << id << ":neighbor_list";

  if (contact_search.getType() == _cst_expli) {
    neighbor_list = new NodesNeighborList(sstr.str());
    nodes_neighbor_list = true;
  }
  else {
    neighbor_list = new NeighborList(sstr.str());
    nodes_neighbor_list = false;
  }

  AKANTU_DEBUG_OUT();
}
