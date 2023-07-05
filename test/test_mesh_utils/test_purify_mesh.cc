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

#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  Mesh mesh(2);

  mesh.read("purify_mesh.msh");

  MeshUtils::purifyMesh(mesh);

  mesh.write("purify_mesh_after.msh");

  if (mesh.getNbNodes() != 21)
    AKANTU_ERROR(
        "The purified mesh does not contain the good number of nodes.");

  if (mesh.getNbElement(_quadrangle_8) != 4)
    AKANTU_ERROR(
        "The purified mesh does not contain the good number of element.");

  return 0;
}
