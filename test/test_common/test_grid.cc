/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 * 
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 * 
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_grid_dynamic.hh"
#include "mesh.hh"
#include "mesh_io.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  const Int spatial_dimension = 2;
  akantu::initialize(argc, argv);

  Mesh circle(spatial_dimension);
  circle.read("circle.msh");

  const auto & l = circle.getLocalLowerBounds();
  const auto & u = circle.getLocalUpperBounds();

  Vector<Real> s{.2, .2};
  Vector<Real> c = (u + l) / 2;

  SpatialGrid<Element> grid(spatial_dimension, s, c);

  Vector<Real> bary(spatial_dimension);
  Element el;
  el.ghost_type = _not_ghost;

  for (const auto & type : circle.elementTypes(spatial_dimension)) {
    Int nb_element = circle.getNbElement(type);
    el.type = type;

    for (Int e = 0; e < nb_element; ++e) {
      el.element = e;
      circle.getBarycenter(el, bary);
      grid.insert(el, bary);
    }
  }

  std::cout << grid << std::endl;
  Mesh mesh(spatial_dimension, "save");
  grid.saveAsMesh(mesh);
  mesh.write("grid.msh");

  akantu::finalize();

  return EXIT_SUCCESS;
}
