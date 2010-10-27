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
#include "mesh.hh"


/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<UInt spatial_dimension> 
class RegularGridNeighborStructure : public ContactNeighborStructure {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RegularGridNeighborStructure(const ContactSearch & contact_search,
			       const Surface & master_surface,
			       const ContactNeighborStructureID & id = "contact_neighbor_structure_id");

  virtual ~RegularGridNeighborStructure();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// initialize the structure
  virtual void init();

  /// update the structure
  virtual void update();

  /// check if an update is needed
  virtual bool check();

  /// function to print the contain of the class
  //virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// compute neighbor cells for a given cell and return number of found neighbor cells
  inline UInt computeNeighborCells(UInt cell, UInt * neighbors, UInt * directional_nb_cells);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set grid spacing
  inline void setGridSpacing(Real spacing, UInt component);

  /// get grid spacing
  inline Real getGridSpacing(UInt component) const;
  

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the mesh
  const Mesh & mesh;

  /// spatial dimension
  //UInt spatial_dimension;

  /// grid spacing
  Real grid_spacing[3];

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
