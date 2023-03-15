/**
 * Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <non_linear_solver.hh>
#include <solid_mechanics_model_cohesive.hh>

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Int spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("bar.msh");

  // Create model
  SolidMechanicsModelCohesive model(mesh);
  model.initFull(_analysis_method = _static, _is_extrinsic = true);

  // Configure solver
  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-10);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  model.initNewSolver(_explicit_lumped_mass);

  // Dynamic insertion of cohesive elements
  model.updateAutomaticInsertion();

  auto dt_crit = model.getStableTimeStep(); // Critical time step
  auto dt = dt_crit * 0.1;                  // Adopted time step
  auto time_simulation = 6.0e-6;            // Total time of simulation (s)
  auto n_steps = int(time_simulation / dt); // Number of time steps

  // Apply Dirichlet BC to block displacements at y direction on top and
  // bottom
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "YBlocked");

  // Constant velocity boundary condition
  // Applied strain rate(s - 1)
  auto strain_rate = 1e5;
  // Applied velocity at the boundary
  auto vel =
      strain_rate * (mesh.getUpperBounds()(_x) - mesh.getLowerBounds()(_x)) / 2;

  model.applyBC(BC::Dirichlet::IncrementValue(-vel * dt, _x), "left");
  model.applyBC(BC::Dirichlet::IncrementValue(vel * dt, _x), "right");

  // VTK plot setup
  model.setBaseName("bar");
  model.addDumpFieldVector("displacement");
  // VTK plot setup for Cohesive model
  model.setBaseNameToDumper("cohesive elements", "cohesive");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");

  for (auto _ [[gnu::unused]] : arange(n_steps)) {
    // Apply velocity at the extremities
    model.applyBC(BC::Dirichlet::IncrementValue(-vel * dt, _x), "left");
    model.applyBC(BC::Dirichlet::IncrementValue(vel * dt, _x), "right");

    model.dump();
    model.dump("cohesive elements");

    // Run simulation
    model.checkCohesiveStress();
    model.solveStep("explicit_lumped");
  }
}
