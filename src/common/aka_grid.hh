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

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void beginInsertions() { data.countToCSR(); data.resizeCols(); data.beginInsertions(); }
  void endInsertions() { data.endInsertions(); }
  void insert(const T & d, const types::RVector & position);

  inline void count(const types::RVector & position);
  inline UInt getNumCell(const types::RVector & position) const;
  inline UInt getCell(Real position, UInt direction) const;

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const {};

  typedef typename CSR<T>::iterator iterator;
  typedef typename CSR<T>::const_iterator const_iterator;
  inline iterator beginCell(UInt cell) { return data.begin(cell); };
  inline iterator endCell(UInt cell) { return data.end(cell); };

  inline const_iterator beginCell(UInt cell) const { return data.begin(cell); };
  inline const_iterator endCell(UInt cell) const { return data.end(cell); };

  struct neighbor_cells_iterator : std::iterator<std::forward_iterator_tag, UInt> {
    neighbor_cells_iterator(const RegularGrid & grid,
			    UInt cell,
			    bool end) : grid(grid) {
      this->cell = cell;
      it = end ? 8 : -1; // if end next cell == 9 else next cell == 0
      nextCell();
    };

    neighbor_cells_iterator(const neighbor_cells_iterator & it) : grid(it.grid) {
      if(this != &it) {
	this->cur_cell = it.cur_cell;
	this->cell = it.cell;
	this->it = it.it;
      }
    };

    neighbor_cells_iterator& operator++() { nextCell(); return *this; };
    neighbor_cells_iterator operator++(int) { neighbor_cells_iterator tmp(*this); operator++(); return tmp; };
    bool operator==(const neighbor_cells_iterator& rhs) { return cell == rhs.cell && it = rhs.it; };
    bool operator!=(const neighbor_cells_iterator& rhs) { return cell != rhs.cell || it != rhs.it; };

    UInt operator*() { return cur_cell; };

  private:
    UInt nextCell() {
      ++it;
      cur_cell = cell;
      UInt t = 1;
      Int position[3] = {0,0,0};
      for (UInt i = 0; i < grid.dimension; ++i, t *= 3) position[i] = ((it / t) % 3) - 1;
      for (UInt i = grid.dimension - 1; i > 0; --i) cur_cell += position[i] * grid.nb_cells[i - 1];
      cur_cell += position[0];
      return cur_cell;
    }

  private:
    UInt cur_cell; //current position of the iterator
    UInt cell; //central cell
    UInt it; // number representing the current neighbor in base 3;

    const RegularGrid & grid;
  };

  inline neighbor_cells_iterator beginNeighborCells(UInt cell) {
    return neighbor_cells_iterator(*this, cell, false);
  }

  inline neighbor_cells_iterator endNeighborCells(UInt cell) {
    return neighbor_cells_iterator(*this, cell, true);
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  inline UInt getNbCells(UInt dim) const { return nb_cells[dim]; };

  inline UInt getCell(UInt num_cells_by_direction[]) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  CSR<T> data;

  UInt dimension;

  Real lower_bounds[3];
  Real upper_bounds[3];

  UInt nb_cells[3];

  // UInt total_nb_cells;

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
