/**
 * @file   aka_grid_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Thu Jul 28 15:27:55 2011
 *
 * @brief  
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

  UInt total_nb_cells = 1;

  std::fill_n(this->upper_bounds, 3, 0);
  std::fill_n(this->lower_bounds, 3, 0);
  std::fill_n(this->spacing, 3, 0);
  std::fill_n(this->nb_cells, 3, 0);

  for (UInt i = 0; i < dimension; ++i) {
    this->upper_bounds[i] = upper_bounds[i];
    this->lower_bounds[i] = lower_bounds[i];
    this->spacing[i]      = spacing[i];

    // +2 to add an extra cell on each side
    this->nb_cells[i] = UInt(std::ceil((upper_bounds[i] - lower_bounds[i]) / spacing[i]) + 2);

    total_nb_cells *= this->nb_cells[i];

    std::cout << nb_cells[i] << std::endl;
    //this->cells_origin[i] = (std::floor(lower_bounds[i] / spacing[i]) - 1) * spacing[i];
  }

  this->data.resizeRows(total_nb_cells);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::insert(const T & d, const types::RVector & position) {
  UInt num_cell = getNumCell(position);
  data.insertInRow(num_cell, d);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::count(const types::RVector & position) {
  UInt num_cell = getNumCell(position);
  data.rowOffset(num_cell)++;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt RegularGrid<T>::getNumCell(const types::RVector & position) const {
  UInt cell[dimension];
  for (UInt i = 0; i < dimension; ++i) {
    cell[i] = getCell(position(i), i);
  }

  UInt num_cell = getCell(cell);

  return num_cell;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt RegularGrid<T>::getCell(Real position, UInt direction) const {
  return UInt(std::floor((position - lower_bounds[direction]) / spacing[direction]) + 1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt RegularGrid<T>::getCell(UInt cell[]) const {
  UInt num_cell = 0;
  for (UInt i = dimension - 1; i > 0; --i) {
    num_cell += cell[i];
    num_cell *= nb_cells[i - 1];
  }
  num_cell += cell[0];

  return num_cell;
}
