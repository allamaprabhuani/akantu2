/**
 * @file   aka_grid_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Aug 31 11:09:48 2011
 *
 * @brief  implementation of template functions of the RegularGrid class
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

/* -------------------------------------------------------------------------- */
template<typename T>
RegularGrid<T>::RegularGrid(UInt dimension, Real * lower_bounds,
			    Real * upper_bounds, Real * spacing) : dimension(dimension) {

  //  UInt total_nb_cells = 1;
  total_nb_cells = 1;
  this->dimension = dimension;

  std::fill_n(this->upper_bounds, 3, 0);
  std::fill_n(this->lower_bounds, 3, 0);
  std::fill_n(this->spacing, 3, 0);
  std::fill_n(this->nb_cells, 3, 0);

  for (UInt i = 0; i < dimension; ++i) {
    // +2 to add an extra cell on each side

    this->nb_cells[i] = UInt(std::ceil((upper_bounds[i] - lower_bounds[i]) / spacing[i])) + 2;
    this->lower_bounds[i] = lower_bounds[i] - spacing[i];
    this->upper_bounds[i] = this->lower_bounds[i] + this->nb_cells[i] * spacing[i];

    this->spacing[i]      = spacing[i];

    total_nb_cells *= this->nb_cells[i];
  }
  this->data.resize(total_nb_cells);
}


/* -------------------------------------------------------------------------- */
template<typename T>
class RegularGrid<T>::Cell {
public:
  Cell() {
    id = 0;
    position[0] = position[1] = position[2] = 0;
    grid = NULL;
  }

  Cell(const RegularGrid<T> & grid) : grid(&grid) { id = 0; std::fill_n(position, 3, 0); }

  Cell(const Cell & cell) {
    if(&cell != this) {
      id = cell.id;
      position[0] = cell.position[0];
      position[1] = cell.position[1];
      position[2] = cell.position[2];
      grid = cell.grid;
    }
  }

  Cell & operator=(const Cell & cell) {
    if(&cell != this) {
      id = cell.id;
      position[0] = cell.position[0];
      position[1] = cell.position[1];
      position[2] = cell.position[2];
      grid = cell.grid;
    }
    return *this;
  }

  bool operator==(const Cell & cell) const { return id == cell.id; }
  bool operator!=(const Cell & cell) const { return id != cell.id; }

  inline void updateID() {
    id = 0;
    for (UInt i = grid->getDimension() - 1; i > 0; --i) {
      id += position[i];
      id *= grid->getNbCells(i - 1);
    }
    id += position[0];
  }

  friend class RegularGrid<T>;
  friend class GridSynchronizer;
private:
  const RegularGrid<T> * grid;

  UInt id;
  UInt position[3];
};

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::iterator RegularGrid<T>::beginCell(const typename RegularGrid<T>::Cell & cell) {
  return data(cell.id).begin();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::iterator RegularGrid<T>::endCell(const typename RegularGrid<T>::Cell cell) {
  return data(cell.id).end();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::const_iterator RegularGrid<T>::beginCell(const typename RegularGrid<T>::Cell & cell) const {
  return data(cell.id).begin();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::const_iterator RegularGrid<T>::endCell(const typename RegularGrid<T>::Cell & cell) const {
  return data(cell.id).end();
}


/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::insert(const T & d, const Vector<Real> & position) {
  Cell cell = getCell(position);
  UInt num_cell = cell.id;
  AKANTU_DEBUG_ASSERT(num_cell < total_nb_cells,
		      "The position of " << d << " is not in the grid (" << num_cell
		      << " > " << total_nb_cells << ") ["
		      << cell.position[0] << ", " << cell.position[1] << ", "<< cell.position[2] << "] : "
		      << position << " -- grid info " << *this);
  data(num_cell).push_back(d);
}

// /* -------------------------------------------------------------------------- */
// template<typename T>
// void RegularGrid<T>::count(const Vector<Real> & position) {
//   Cell cell = getCell(position);
//   UInt num_cell = cell.id;
//   // std::cout << num_cell << " - "
//   // 	    << cell.position[0] << ", " << cell.position[1] << ", " << cell.position[2] << " : "
//   // 	    << position[0] << ", " << position[1] << ", " << position[2]
//   // 	    << std::endl;
//   data.rowOffset(num_cell)++;
// }

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::Cell RegularGrid<T>::getCell(const Vector<Real> & position) const {
  Cell cell(*this);
  for (UInt i = 0; i < dimension; ++i) {
    cell.position[i] = getCell(position(i), i);
  }

  cell.updateID();

  return cell;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt RegularGrid<T>::getCell(Real position, UInt direction) const {
  return UInt(std::floor((position - lower_bounds[direction]) / spacing[direction]));
}

/* -------------------------------------------------------------------------- */
template<typename T>
struct RegularGrid<T>::neighbor_cells_iterator : private std::iterator<std::forward_iterator_tag, UInt> {
  neighbor_cells_iterator(const RegularGrid & grid,
			  const Cell & cell,
			  bool end) : cell(cell), grid(grid) {
    std::fill_n(position_start, 3, -1);
    std::fill_n(position_end  , 3,  1);

    for (UInt i = 0; i < 3; ++i) {
      if((grid.getNbCells(i) == 0) ||                  // no cells in this direction
	 (cell.position[i] == 0))                      // first cell in this direction
	position_start[i] = 0;

      if((grid.getNbCells(i) == 0) ||                  // no cells in this direction
	 (cell.position[i] == grid.getNbCells(i) - 1)) // last cell in this direction
	position_end[i] = 0;
      position[i] = end ? position_end[i] : position_start[i];
    }

    updateIt();
    if(end) { (this->it)++; }
  };

  neighbor_cells_iterator(const neighbor_cells_iterator & it) :  cell(it.cell),
								 it(it.it),
								 grid(it.grid) {
    std::copy(it.position      , it.position       + 3, position      );
    std::copy(it.position_start, it.position_start + 3, position_start);
    std::copy(it.position_end  , it.position_end   + 3, position_end  );
  };

  neighbor_cells_iterator& operator++() {
    bool last = false;
    if(position[0] < position_end[0]) {
      position[0]++;
    } else {
      position[0] = position_start[0];
      if(position[1] < position_end[1]) {
	position[1]++;
      } else {
	position[1] = position_start[1];
	if(position[2] < position_end[2]) {
	  position[2]++;
	} else {
	  last = true;
	}
      }
    }

    if(last) ++it;
    else {
      updateIt();
    }
    return *this;
  };

  neighbor_cells_iterator operator++(int) { neighbor_cells_iterator tmp(*this); operator++(); return tmp; };
  bool operator==(const neighbor_cells_iterator& rhs) { return cell == rhs.cell && it = rhs.it; };
  bool operator!=(const neighbor_cells_iterator& rhs) { return cell != rhs.cell || it != rhs.it; };

  Cell operator*() {
    Cell cur_cell(grid);
    for (UInt i = 0; i < 3; ++i) cur_cell.position[i] = cell.position[i] + position[i];
    cur_cell.updateID();
    return cur_cell;
  };

private:
  void updateIt() {
    it = 0;
    for (UInt i = 0; i < 3; ++i) it = it * 3 + (position[i] + 1);
  }

private:
  //  UInt cur_cell; //current position of the iterator
  const Cell & cell; //central cell
  UInt it; // number representing the current neighbor in base 3;

  Int position_start[3], position_end[3], position[3];

  const RegularGrid & grid;
};

/* -------------------------------------------------------------------------- */
template<typename T>
inline typename RegularGrid<T>::neighbor_cells_iterator RegularGrid<T>::beginNeighborCells(const typename RegularGrid<T>::Cell & cell) {
  return neighbor_cells_iterator(*this, cell, false);
}

template<typename T>
inline typename RegularGrid<T>::neighbor_cells_iterator RegularGrid<T>::endNeighborCells(const typename RegularGrid<T>::Cell & cell) {
  return neighbor_cells_iterator(*this, cell, true);
}


/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::printself(std::ostream & stream, int indent) const {
  std::string space;
  for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  std::streamsize prec        = stream.precision();
  std::ios_base::fmtflags ff  = stream.flags();

  stream.setf (std::ios_base::showbase);
  stream.precision(5);


  stream << space << "RegularGrid<" << debug::demangle(typeid(T).name()) << "> [" << std::endl;
  stream << space << " + dimension    : " << this->dimension << std::endl;
  stream << space << " + lower_bounds : {"
	 << this->lower_bounds[0] << ", "
	 << this->lower_bounds[1] << ", "
	 << this->lower_bounds[2] << "}" << std::endl;
  stream << space << " + upper_bounds : {"
	 << this->upper_bounds[0] << ", "
	 << this->upper_bounds[1] << ", "
	 << this->upper_bounds[2] << "}" << std::endl;
  stream << space << " + spacing      : {"
	 << this->spacing[0] << ", "
	 << this->spacing[1] << ", "
	 << this->spacing[2] << "}" << std::endl;
  stream << space << " + nb_cells     : " << this->total_nb_cells << " - {"
	 << this->nb_cells[0] << ", "
	 << this->nb_cells[1] << ", "
	 << this->nb_cells[2] << "}" << std::endl;
  stream << space << "]" << std::endl;

  stream.precision(prec);
  stream.flags(ff);

}


/* -------------------------------------------------------------------------- */
__END_AKANTU__

#include "mesh.hh"

__BEGIN_AKANTU__

template<typename T>
void RegularGrid<T>::saveAsMesh(Mesh & mesh) const {
  Array<Real> & nodes = const_cast<Array<Real> &>(mesh.getNodes());
  UInt nb_nodes = 1;


  for (UInt i = 0; i < dimension; ++i) {
    nb_nodes *= (nb_cells[i] + 1);
  }
  nodes.resize(nb_nodes);

  if(dimension == 1) {
    for (UInt n = 0; n < nb_nodes; ++n) {
      nodes(n, 0) = n * spacing[0] + lower_bounds[0];
    }

    mesh.addConnectivityType(_segment_2);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(_segment_2));
    connectivity.resize(total_nb_cells);

    for (UInt e = 0; e < nb_cells[0]; ++e) {
      connectivity(e, 0) = e;
      connectivity(e, 1) = e + 1;
    }
  }

  if(dimension == 2) {
    UInt nnx = nb_cells[0] + 1;
    UInt nny = nb_cells[1] + 1;

    for (UInt nx = 0; nx < nnx; ++nx) {
      for (UInt ny = 0; ny < nny; ++ny) {
	UInt n = nx * nny + ny;
	nodes(n, 0) = nx * spacing[0] + lower_bounds[0];
	nodes(n, 1) = ny * spacing[1] + lower_bounds[1];
      }
    }

    mesh.addConnectivityType(_quadrangle_4);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(_quadrangle_4));
    connectivity.resize(total_nb_cells);
    for (UInt ex = 0; ex < nb_cells[0]; ++ex) {
      for (UInt ey = 0; ey < nb_cells[1]; ++ey) {
	UInt e = (ex * nb_cells[1] + ey);
	connectivity(e, 0) = ex       * nny + ey;
	connectivity(e, 1) = (ex + 1) * nny + ey;
	connectivity(e, 2) = (ex + 1) * nny + ey + 1;
	connectivity(e, 3) = ex       * nny + ey + 1;
      }
    }
  }

  if(dimension == 3) {
    UInt nnx = nb_cells[0] + 1;
    UInt nny = nb_cells[1] + 1;
    UInt nnz = nb_cells[2] + 1;

    for (UInt nx = 0; nx < nnx; ++nx) {
      for (UInt ny = 0; ny < nny; ++ny) {
	for (UInt nz = 0; nz < nnz; ++nz) {
	  UInt n = (nx * nny + ny) * nnz + nz;
	  nodes(n, 0) = nx * spacing[0] + lower_bounds[0];
	  nodes(n, 1) = ny * spacing[1] + lower_bounds[1];
	  nodes(n, 2) = nz * spacing[2] + lower_bounds[2];
	}
      }
    }

    mesh.addConnectivityType(_hexahedron_8);
    Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh.getConnectivity(_hexahedron_8));
    connectivity.resize(total_nb_cells);
    for (UInt ex = 0; ex < nb_cells[0]; ++ex) {
      for (UInt ey = 0; ey < nb_cells[1]; ++ey) {
	for (UInt ez = 0; ez < nb_cells[2]; ++ez) {
	  UInt e = (ex * nb_cells[1] + ey) * nb_cells[2] + ez;
	  connectivity(e, 0) = (ex     * nny + ey  ) * nnz + ez;
	  connectivity(e, 1) = ((ex+1) * nny + ey  ) * nnz + ez;
	  connectivity(e, 2) = ((ex+1) * nny + ey+1) * nnz + ez;
	  connectivity(e, 3) = (ex     * nny + ey+1) * nnz + ez;
	  connectivity(e, 4) = (ex     * nny + ey  ) * nnz + ez+1;
	  connectivity(e, 5) = ((ex+1) * nny + ey  ) * nnz + ez+1;
	  connectivity(e, 6) = ((ex+1) * nny + ey+1) * nnz + ez+1;
	  connectivity(e, 7) = (ex     * nny + ey+1) * nnz + ez+1;
	}
      }
    }
  }
}
