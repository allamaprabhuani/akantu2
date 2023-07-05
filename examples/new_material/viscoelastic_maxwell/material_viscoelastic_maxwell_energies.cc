/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
/* -------------------------------------------------------------------------- */
#include "material_viscoelastic_maxwell.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"

using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {
  akantu::initialize("material_viscoelastic_maxwell.dat", argc, argv);

  // sim data
  Real eps = 0.1;

  const Int dim = 2;
  Real sim_time = 100.;
  Real T = 10.;
  Mesh mesh(dim);
  mesh.read("material_viscoelastic_maxwell_mesh.msh");

  SolidMechanicsModel model(mesh);

  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initFull(_analysis_method = _static);
  std::cout << model.getMaterial(0) << std::endl;

  std::stringstream filename_sstr;
  filename_sstr << "material_viscoelastic_maxwell_output.out";
  std::ofstream output_data;
  output_data.open(filename_sstr.str().c_str());

  Material & mat = model.getMaterial(0);

  Real time_step = 0.1;

  const Array<Real> & coordinate = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & blocked = model.getBlockedDOFs();

  /// Setting time step

  model.setTimeStep(time_step);

  model.setBaseName("dynamic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("grad_u");
  model.addDumpField("stress");
  model.addDumpField("strain");

  Int max_steps = sim_time / time_step + 1;
  Real time = 0.;

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 10);
  solver.set("threshold", 1e-7);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for (Int s = 0; s <= max_steps; ++s) {

    std::cout << "Time Step = " << time_step << "s" << std::endl;
    std::cout << "Time = " << time << std::endl;

    // impose displacement
    Real epsilon = 0;
    if (time < T) {
      epsilon = eps * time / T;
    } else {
      epsilon = eps;
    }

    for (auto && [coord, disp, block] :
         zip(make_view(coordinate, dim), make_view(displacement, dim),
             make_view(blocked, dim))) {
      if (Math::are_float_equal(coord(_x), 0.0)) {
        disp(_x) = 0;
        disp(_y) = epsilon * coord(_y);
        block.set(true);
      } else if (Math::are_float_equal(coord(_y), 0.0)) {
        disp(_x) = epsilon * coord(_x);
        disp(_y) = 0;
        block.set(true);
      } else if (Math::are_float_equal(coord(_x), 0.001)) {
        disp = epsilon * coord;
        block.set(true);
      } else if (Math::are_float_equal(coord(_y), 0.001)) {
        disp = epsilon * coord;
        block.set(true);
      }
    }

    try {
      model.solveStep();
    } catch (debug::NLSNotConvergedException & e) {
      std::cout << "Didn't converge after " << e.niter
                << " iterations. Error is " << e.error << std::endl;
      return EXIT_FAILURE;
    }

    // for debugging
    // auto int_force = model.getInternalForce();
    // auto &K = model.getDOFManager().getMatrix("K");
    // K.saveMatrix("K.mtx");

    Int nb_iter = solver.get("nb_iterations");
    std::cout << "Converged in " << nb_iter << " iterations" << std::endl;

    model.dump();

    Real epot = mat.getEnergy("potential");
    Real edis = mat.getEnergy("dissipated");
    Real work = mat.getEnergy("work");

    // data output
    output_data << s * time_step << " " << epsilon << " " << epot << " " << edis
                << " " << work << std::endl;
    time += time_step;
  }
  output_data.close();
}
