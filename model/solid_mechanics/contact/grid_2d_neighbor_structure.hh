/**
 * @file   grid_2d_neighbor_structure.hh
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @date   Tue Dec  7 12:54:48 2010
 *
 * @brief  Class which creates the neighbor lists (with a grid) to handle contact in 2d 
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

#ifndef __AKANTU_GRID_2D_NEIGHBOR_STRUCTURE_HH__
#define __AKANTU_GRID_2D_NEIGHBOR_STRUCTURE_HH__

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
class Grid2dNeighborStructure : public ContactNeighborStructure{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Grid2dNeighborStructure(const ContactSearch & contact_search,
			  const Surface & master_surface,
			  const ContactNeighborStructureType & type,
			  const ContactNeighborStructureID & id = "contact_neighbor_structure_id");
  virtual ~Grid2dNeighborStructure();
  
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
  /// main routine to create the grid and fill the neighbor list
  void createGrid(bool initial_position);

  /// get extrema of each surface
  void getBounds(Real * coord, Real * x_bounds, Real * y_bounds);

  /// get intersection space between master surface and slave ones
  bool getBoundsIntersection(Real *x_bounds, Real * y_bounds, Real * x_int, Real * y_int);

  /// get minimum mesh size (segment length)
  Real getMinSize(Real * coord);

  /// which grid cells are intersected by segment of the master surface
  void traceSegments(Real * coord, Real * origin, UInt * nb_cells,
		     UInt * cell_to_seg_off, Vector<UInt> & cell_to_segments);

  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Mesh
  Mesh & mesh;

  /// grid spacing
  Real spacing;

  /// maximal displacement for grid update
  Real max_increment[2];
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "2d_grid_neighbor_structure_inline_impl.cc"

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const Grid2dNeighborStructure & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__



#endif /* __AKANTU_GRID_2D_NEIGHBOR_STRUCTURE_HH__ */
