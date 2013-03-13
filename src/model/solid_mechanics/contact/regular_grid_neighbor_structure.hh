/**
 * @file   regular_grid_neighbor_structure.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Mon Oct 11 18:09:29 2010
 *
 * @brief  Structure that handles the neighbor lists by a regular grid
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

#ifndef __AKANTU_REGULAR_GRID_NEIGHBOR_STRUCTURE_HH__
#define __AKANTU_REGULAR_GRID_NEIGHBOR_STRUCTURE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_vector.hh"
#include "contact_neighbor_structure.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class NodesNeighborList : public NeighborList {
public:
  NodesNeighborList(const ID & id);
  virtual ~NodesNeighborList() {};
public:
  ID id;
  /// neighbor master nodes
  Array<UInt> master_nodes_offset;
  Array<UInt> master_nodes;
};


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class RegularGridNeighborStructure : public ContactNeighborStructure {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RegularGridNeighborStructure(const ContactSearch & contact_search,
			       const Surface & master_surface,
			       const ContactNeighborStructureType & type,
			       const ContactNeighborStructureID & id = "contact_neighbor_structure_id");

  virtual ~RegularGridNeighborStructure();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the structure
  void initNeighborStructure();

  /// update the structure
  void update();

  /// check if an update is needed
  bool check();

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// compute neighbor structure
  void update(Real * node_position);

  /// construct neighbor list
  void constructNeighborList(Int directional_nb_cells[spatial_dimension],
			     UInt nb_cells,
			     Array<Int> * cell,
			     UInt * impactor_nodes_cell_offset,
			     UInt * impactor_nodes_cell,
			     UInt * master_nodes_cell_offset,
			     UInt * master_nodes_cell);


  /// construct nodes neighbor list
  void constructNodesNeighborList(Int directional_nb_cells[spatial_dimension],
				  UInt nb_cells,
				  Array<Int> * cell,
				  UInt * impactor_nodes_cell_offset,
				  UInt * impactor_nodes_cell,
				  UInt * master_nodes_cell_offset,
				  UInt * master_nodes_cell);


  /// compute neighbor cells for a given cell and return number of found neighbor cells
  inline UInt computeNeighborCells(UInt cell, UInt * neighbors, Int * directional_nb_cells);

  /// compute global cell number given the directional cell number
  inline UInt computeCellNb(Int * directional_nb_cells, Int * directional_cell);

  /// initializes the neighbor list
  inline void constructNeighborList();

  /// compute minimal grid size and set it
  void setMinimalGridSpacing();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set grid spacing
  inline void setGridSpacing(Real spacing, UInt component);

  /// get grid spacing
  inline Real getGridSpacing(UInt component) const;

  /// set security factor
  inline void setSecurityFactor(Real factor, UInt component);

  /// get security factor
  inline Real getSecurityFactor(UInt component) const;

  /// set max increment
  inline void setMaxIncrement(Real increment, UInt component);

  /// get max increment
  inline Real getMaxIncrement(UInt component) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the mesh
  const Mesh & mesh;

  /// type of neighbor list to create
  bool nodes_neighbor_list;

  /// grid spacing
  Real grid_spacing[3];

  /// maximal displacement since last grid update
  Real max_increment[3];

  /// security factor for grid update test
  Real security_factor[3];
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "regular_grid_neighbor_structure_inline_impl.cc"
#endif

/// standard output stream operator
/*inline std::ostream & operator <<(std::ostream & stream, const RegularGridNeighborStructure & _this)
{
  _this.printself(stream);
  return stream;
  }*/


__END_AKANTU__

#endif /* __AKANTU_REGULAR_GRID_NEIGHBOR_STRUCTURE_HH__ */
