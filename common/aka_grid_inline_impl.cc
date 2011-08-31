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
    this->nb_cells[i]     = std::ceil((upper_bounds[i] - lower_bounds[i]) / spacing[i]) + 2;

    total_nb_cells *= this->nb_cells[i];
    //this->cells_origin[i] = (std::floor(lower_bounds[i] / spacing[i]) - 1) * spacing[i];
  }

  this->data.resizeRows(total_nb_cells);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::insert(const T & d, const types::RVector & position) {
  UInt num_cell = getCell(position);
  data.insertInRow(num_cell, d);
}

/* -------------------------------------------------------------------------- */
template<typename T>
void RegularGrid<T>::count(const types::RVector & position) {
  UInt num_cell = getCell(position);
  data.rowOffset(num_cell)++;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline UInt RegularGrid<T>::getCell(const types::RVector & position) const {
  UInt cell[dimension];
  for (UInt i = 0; i < dimension; ++i) {
    cell[i] = std::floor(position(i) / spacing[i]) + 1;
  }

  UInt num_cell = 0;
  for (UInt i = dimension - 1; i > 0; --i) {
    num_cell += cell[i];
    num_cell *= nb_cells[i - 1];
  }
  num_cell += cell[0];

  return num_cell;
}

// /* -------------------------------------------------------------------------- */
// __END_AKANTU__
// #include "mesh.hh"
// #include "mesh_io.hh"
// __BEGIN_AKANTU__

// template<typename T>
// inline UInt RegularGrid<T>::saveAsMesh() const {
//   Mesh mesh(spatial_dimension);

//   Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());

//   ElementType type;
//   switch(spatial_dimension) {
//   case 1: type = _segment_2; break;
//   case 2: type = _quadrangle_4; break;
//   case 3: type = _hexahedron_8; break;
//   }

//   mesh.addConnecticityType(type);
//   Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(mesh.getConnectivity(type));

//   connectivity.resize(data.getNbRows());

//   UInt nb_nodes = 1;
//   for (UInt i = 0; i < spatial_dimension; ++i) nb_nodes *= nb_cells[i] + 1;
//   nodes.resize(nb_nodes);

//   UInt n = 0;
//   for (UInt i = 0; i < nb_cells[0]; ++i) {
//     nodes(n, 0) = lower_bounds[0] + (i - 1) * spacing[0];
//     if(spatial_dimension >= 2) {
//       for (UInt j = 0; j < nb_cells[1]; ++j) {
// 	nodes(n, 1) = lower_bounds[1] + (j - 1) * spacing[1];
// 	if(spatial_dimension >= 3) {
// 	  for (UInt k = 0; k < nb_cells[2]; ++k, ++n)
// 	    nodes(n, 2) = lower_bounds[2] + (k - 1) * spacing[2];
// 	} else ++n;
//       }
//     } else ++n;
//   }
//  /// need the connectivity
// }
