/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "heat_transfer_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
Int spatial_dimension = 2;
std::string base_name;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square.msh");

  HeatTransferModel model(mesh);
  // initialize everything
  model.initFull(_analysis_method = _static);

  // boundary conditions
  const Array<Real> & nodes = mesh.getNodes();
  Array<bool> & blocked_dofs = model.getBlockedDOFs();
  Array<Real> & temperature = model.getTemperature();
  double length = 1.;
  Int nb_nodes = nodes.size();
  for (Int i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    Real dx = nodes(i, 0);
    Real dy = nodes(i, 1);

    Vector<Real> dX = {dx, dy};
    dX.array() -= length / 4.;
    Real d = dX.norm();
    if (d < 0.1) {
      blocked_dofs(i) = true;
      temperature(i) = 300.;
    }

    if (std::abs(dx) < 1e-4 || std::abs(dy) < 1e-4)
      blocked_dofs(i) = true;
    if (std::abs(dx - length) < 1e-4 || std::abs(dy - length) < 1e-4)
      blocked_dofs(i) = true;
  }

  model.setBaseName("heat_diffusion_static_2d");
  model.addDumpField("temperature");
  model.addDumpField("internal_heat_rate");
  model.addDumpField("conductivity");
  model.addDumpField("blocked_dofs");
  model.dump();

  model.solveStep();
  model.dump();

  return 0;
}
