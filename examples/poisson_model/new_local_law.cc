/**
 * @file   diffusion_dynamics_2d.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Fri Mar 16 2018
 *
 * @brief  Example of diffusion constitutive law
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "poisson_model.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 1;
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  initialize("local_law.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("bar.msh");

  PoissonModel model(mesh);

  // initialize everything
  model.initFull();

  // get stable time step
  Real time_step = model.getStableTimeStep() * 0.1;
  std::cout << "time step is:" << time_step << std::endl;
  model.setTimeStep(time_step);

  // boundary conditions
  const Array<Real> & nodes = model.getFEEngine().getMesh().getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & concentration = model.getDof();

  auto & external_flux_rate = model.getExternalDofRate();
  external_flux_rate(0) = 1e-8;

  std::cout << external_flux_rate;

  model.setBaseName("new_local1d");
  model.addDumpField("dof");
  model.addDumpField("internal_dof_rate");
  model.addDumpField("external_dof_rate");

  model.dump();
  
  // main loop
  int max_steps = 15000;
  for (int i = 0; i < max_steps; i++) {
    model.solveStep();

    if (i % 100 == 0)
      model.dump();

    std::cout << "Step " << i << "/" << max_steps << std::endl;
  }
  std::cout << "\n\n Stable Time Step is : " << time_step << "\n \n"
            << std::endl;


  return 0;
}
