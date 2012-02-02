/**
 * @file   aka_grid.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 28 11:44:18 2011
 *
 * @brief  A regular grid that can contain information per cells
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
#include "aka_common.hh"
#include "aka_csr.hh"

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_GRID_HH__
#define __AKANTU_AKA_GRID_HH__

__BEGIN_AKANTU__

class Mesh;

template<typename T>
class RegularGrid {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  RegularGrid() {};

  RegularGrid(UInt dimension, Real * lower_bounds,
	      Real * upper_bounds, Real * spacing);
  virtual ~RegularGrid() {};

  class Cell;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void saveAsMesh(Mesh & mesh) const;

  void beginInsertions() { data.countToCSR(); data.resizeCols(); data.beginInsertions(); }
  void endInsertions() { data.endInsertions(); }
  void insert(const T & d, const types::RVector & position);

  inline void count(const types::RVector & position);
  inline Cell getCell(const types::RVector & position) const;
  inline UInt getCell(Real position, UInt direction) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  typedef typename CSR<T>::iterator iterator;
  typedef typename CSR<T>::const_iterator const_iterator;

  inline iterator beginCell(const Cell & cell);
  inline iterator endCell(const Cell cell);
  inline const_iterator beginCell(const Cell & cell) const;
  inline const_iterator endCell(const Cell & cell) const;

  class neighbor_cells_iterator;

  inline neighbor_cells_iterator beginNeighborCells(const Cell & cell);
  inline neighbor_cells_iterator endNeighborCells(const Cell & cell);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  inline UInt getNbCells(UInt dim) const { return nb_cells[dim]; };

  AKANTU_GET_MACRO(Dimension, dimension, UInt)

  friend class Cell;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  CSR<T> data;

  UInt dimension;

  Real lower_bounds[3];
  Real upper_bounds[3];

  UInt nb_cells[3];

  UInt total_nb_cells;

  Real spacing[3];
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "aka_grid_tmpl.hh"


/// standard output stream operator
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const RegularGrid<T> & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_AKA_GRID_HH__ */
