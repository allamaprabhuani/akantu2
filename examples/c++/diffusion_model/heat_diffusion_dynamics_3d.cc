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
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */

const Int spatial_dimension = 3;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("cube.msh");

  HeatTransferModel model(mesh);
  // initialize everything
  model.initFull(_analysis_method = _explicit_lumped_mass);

  // get and set stable time step
  Real time_step = model.getStableTimeStep() * 0.8;
  std::cout << "Stable Time Step is : " << time_step / .8 << "\n";
  std::cout << "time step is:" << time_step << "\n";
  model.setTimeStep(time_step);

  /// boundary conditions
  const Array<Real> & nodes = mesh.getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & temperature = model.getTemperature();
  auto nb_nodes = mesh.getNbNodes();

  double length = 1.;

  for (Int i = 0; i < nb_nodes; ++i) {
    temperature(i) = 100.;

    // to insert a heat source
    Real dx = nodes(i, 0) - length / 2.;
    Real dy = nodes(i, 1) - length / 2.;
    Real dz = nodes(i, 2) - length / 2.;

    Real d = sqrt(dx * dx + dy * dy + dz * dz);

    if (d < 0.1) {
      boundary(i) = true;
      temperature(i) = 300.;
    }
  }

  model.setBaseName("heat_diffusion_cube3d");
  model.addDumpField("temperature");
  model.addDumpField("temperature_rate");
  model.addDumpField("internal_heat_rate");

  // //for testing
  int max_steps = 1000;

  for (int i = 0; i < max_steps; i++) {
    model.solveStep();

    if (i % 100 == 0) {
      model.dump();
    }

    if (i % 10 == 0) {
      std::cout << "Step " << i << "/" << max_steps << "\n";
    }
  }

  return 0;
}
