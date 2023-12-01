/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "element_group.hh"
#include "mesh_iterators.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const Int spatial_dimension = 2;
  const Int max_steps = 1;

  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  SolidMechanicsModelCohesive model(mesh);
  model.getElementInserter().setLimit(_x, -0.05, 0.05);

  /// model initialization
  model.initFull(_analysis_method = _static
                 ,_is_extrinsic = false
                 );

  auto & solver = model.getNonLinearSolver();
//  solver.set("convergence_type", SolveConvergenceCriteria::_residual);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  solver.set("max_iterations", 2);

  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "point");

  model.setBaseName("cohesive");
  model.addDumpFieldVector("displacement");
//  model.addDumpField("velocity");
//  model.addDumpField("acceleration");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.dump();


  Real increment = 0.01;

  /// Main loop
  for (Int s = 1; s <= max_steps; ++s) {
    model.applyBC(BC::Dirichlet::IncrementValue(increment, _x), "right");
    // debug::setDebugLevel(dblDump);
    debug::setDebugLevel(dblInfo);
    model.solveStep();
    debug::setDebugLevel(dblWarning);

    if (s % 1 == 0) {
      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }
  }

  Real E = model.getEnergy("potential");

  std::cout << "Elastic energy : " << E << std::endl;

//  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
//    std::cout << "The dissipated energy is incorrect" << std::endl;
//    return EXIT_FAILURE;
//  }

  finalize();

  return EXIT_SUCCESS;
}

