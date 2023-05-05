/**
 * Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "coupler_solid_phasefield.hh"
#include "material.hh"
#include "material_phasefield.hh"
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel &, Real &);
/* -------------------------------------------------------------------------- */

int main(int argc, char * argv[]) {

  std::ofstream os("data.csv");
  os << "#strain stress damage analytical_sigma analytical_damage" << std::endl;

  initialize("material_hybrid.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("test_one_element.msh");

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  model.initFull(_analysis_method = _static);
  auto & solver = model.getNonLinearSolver("static");
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-6);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);
  phase.initFull(_analysis_method = _static);
  auto & solver_phase = phase.getNonLinearSolver("static");
  solver_phase.set("max_iterations", 1000);
  solver_phase.set("threshold", 1e-6);
  solver_phase.set("convergence_type", SolveConvergenceCriteria::_residual);

  model.setBaseName("phase_solid");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  UInt nbSteps = 1500;
  Real increment = 1e-4;

  auto & stress = model.getMaterial(0).getArray<Real>("stress", _quadrangle_4);
  auto & damage = model.getMaterial(0).getArray<Real>("damage", _quadrangle_4);

  Real analytical_damage{0.};
  Real analytical_sigma{0.};

  auto & phasefield = phase.getPhaseField(0);

  const Real E = phasefield.getParam("E");
  const Real nu = phasefield.getParam("nu");
  Real c22 = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

  const Real gc = phasefield.getParam("gc");
  const Real l0 = phasefield.getParam("l0");

  Real error_stress{0.};

  Real error_damage{0.};

  Real max_strain{0.};

  for (UInt s = 0; s < nbSteps; ++s) {
    Real axial_strain{0.};
    if (s < 500) {
      axial_strain = increment * s;
    } else if (s < 1000) {
      axial_strain = (1500 - 2 * double(s)) * increment;
    } else {
      axial_strain = (3 * double(s) - 3500) * increment;
    }
    applyDisplacement(model, axial_strain);

    if (axial_strain > max_strain) {
      max_strain = axial_strain;
    }

    coupler.solve("static", "static");
    phase.savePreviousDamage();

    analytical_damage = max_strain * max_strain * c22 /
                        (gc / l0 + max_strain * max_strain * c22);
    if (axial_strain < 0.) {
      analytical_sigma = c22 * axial_strain;
    } else {
      analytical_sigma = c22 * axial_strain * (1 - analytical_damage) *
                         (1 - analytical_damage);
    }

    error_stress =
        std::abs(analytical_sigma - stress(0, 3)) / std::abs(analytical_sigma);

    error_damage = std::abs(analytical_damage - damage(0)) / analytical_damage;

    if ((error_damage > 1e-8 or error_stress > 1e-8) and
        std::abs(axial_strain) > 1e-13) {
      std::cerr << std::left << std::setw(15) << "Step: " << s << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Axial strain: " << axial_strain << std::endl;
      std::cerr << std::left << std::setw(15)
                << "An. damage: " << analytical_damage << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Damage: " << damage(0) << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error damage: " << error_damage << std::endl;
      std::cerr << std::left << std::setw(15)
                << "Error stress: " << error_stress << std::endl;
      return EXIT_FAILURE;
    }

    os << axial_strain << " " << stress(0, 3) << " " << damage(0) << " "
       << analytical_sigma << " " << analytical_damage << " " << error_stress
       << " " << error_damage << std::endl;

    model.dump();
  }

  os.close();
  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void applyDisplacement(SolidMechanicsModel & model, Real & increment) {
  auto & displacement = model.getDisplacement();

  auto & positions = model.getMesh().getNodes();
  auto & blocked_dofs = model.getBlockedDOFs();

  for (Idx n = 0; n < model.getMesh().getNbNodes(); ++n) {
    if (positions(n, 1) == -0.5) {
      displacement(n, 0) = 0;
      displacement(n, 1) = 0;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    } else {
      displacement(n, 0) = 0;
      displacement(n, 1) = increment;
      blocked_dofs(n, 0) = true;
      blocked_dofs(n, 1) = true;
    }
  }
}

/* -------------------------------------------------------------------------- */
