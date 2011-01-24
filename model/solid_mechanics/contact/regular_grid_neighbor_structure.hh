/**
 * @file   regular_grid_neighbor_structure.hh
 * @author David Kammer <david.kammer@epfl.ch>
 * @date   Mon Oct 11 10:35:04 2010
 *
 * @brief  Structure that handles the neighbor lists by a regular grid 
 *
 * @section LICENSE
 *
 * \<insert license here\>
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
namespace akantu {
  class mesh;
}

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
class NodesNeighborList : public NeighborList {
public:
  NodesNeighborList();
  virtual ~NodesNeighborList() {};
public:
  /// neighbor master nodes
  Vector<UInt> master_nodes_offset;
  Vector<UInt> master_nodes;
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
  void constructNeighborList(UInt directional_nb_cells[spatial_dimension], 
			     UInt nb_cells, 
			     Vector<Int> * cell, 
			     UInt * impactor_nodes_cell_offset, 
			     UInt * impactor_nodes_cell, 
			     UInt * master_nodes_cell_offset, 
			     UInt * master_nodes_cell);


  /// construct nodes neighbor list 
  void constructNodesNeighborList(UInt directional_nb_cells[spatial_dimension], 
				  UInt nb_cells, 
				  Vector<Int> * cell, 
				  UInt * impactor_nodes_cell_offset, 
				  UInt * impactor_nodes_cell, 
				  UInt * master_nodes_cell_offset, 
				  UInt * master_nodes_cell);


  /// compute neighbor cells for a given cell and return number of found neighbor cells
  inline UInt computeNeighborCells(UInt cell, UInt * neighbors, UInt * directional_nb_cells);

  /// compute global cell number given the directional cell number
  inline UInt computeCellNb(UInt * directional_nb_cells, Int * directional_cell);

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

#include "regular_grid_neighbor_structure_inline_impl.cc"

/// standard output stream operator
/*inline std::ostream & operator <<(std::ostream & stream, const RegularGridNeighborStructure & _this)
{
  _this.printself(stream);
  return stream;
  }*/


__END_AKANTU__

#endif /* __AKANTU_REGULAR_GRID_NEIGHBOR_STRUCTURE_HH__ */
