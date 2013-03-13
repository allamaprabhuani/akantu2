/**
 * @file   test_grid.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Aug 06 14:09:13 2012
 *
 * @brief  Test the grid object
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "mesh.hh"
#include "mesh_io.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  const UInt spatial_dimension = 2;
  akantu::initialize(argc, argv);

  Mesh circle(spatial_dimension);
  circle.read("circle.msh");

  Real lower[spatial_dimension];
  Real upper[spatial_dimension];

  circle.computeBoundingBox();
  circle.getLocalLowerBounds(lower);
  circle.getLocalUpperBounds(upper);

  Real spacing[spatial_dimension] = {0.2, 0.2};

  Vector<Real> l(lower, spatial_dimension);
  Vector<Real> u(upper, spatial_dimension);

  Vector<Real> s(spacing, spatial_dimension);

  Vector<Real> c = u;
  c += l;
  c /= 2.;

  SpatialGrid<Element> grid(spatial_dimension, s, c);

  Vector<Real> bary(spatial_dimension);
  Element el;
  el.ghost_type = _not_ghost;

  Mesh::type_iterator it        = circle.firstType(spatial_dimension);
  Mesh::type_iterator last_type = circle.lastType (spatial_dimension);
  for(; it != last_type; ++it) {
    UInt nb_element = circle.getNbElement(*it);
    el.type = *it;

    for (UInt e = 0; e < nb_element; ++e) {
      circle.getBarycenter(e, el.type, bary.storage());
      el.element = e;
      grid.insert(el, bary);
    }
  }

  std::cout << grid << std::endl;
  Mesh mesh(spatial_dimension, "save");
  grid.saveAsMesh(mesh);
  mesh.write("grid.msh");

  Vector<Real> pos(spatial_dimension);

  const SpatialGrid<Element>::CellID & id = grid.getCellID(pos);

  SpatialGrid<Element>::neighbor_cells_iterator nit = grid.beginNeighborCells(id);
  SpatialGrid<Element>::neighbor_cells_iterator nend = grid.endNeighborCells(id);

  for(;nit != nend; ++nit) {
    std::cout << std::endl;
    const SpatialGrid<Element>::Cell & cell = grid.getCell(*nit);
    SpatialGrid<Element>::Cell::const_iterator cit = cell.begin();
    SpatialGrid<Element>::Cell::position_iterator pit = cell.begin_pos();
    SpatialGrid<Element>::Cell::const_iterator cend = cell.end();
    for (; cit != cend; ++cit, ++pit) {
      std::cout << *cit << " " << *pit << std::endl;
    }
  }



  akantu::finalize();

  return EXIT_SUCCESS;
}
