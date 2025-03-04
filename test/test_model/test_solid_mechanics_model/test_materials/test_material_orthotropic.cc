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
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  //  akantu::initialize("orthotropic.dat", argc, argv);
  akantu::initialize("orthotropic.dat", argc, argv);
  Int max_steps = 1000;
  Real epot, ekin;

  Mesh mesh(2);
  mesh.read("square.msh");
  mesh.createBoundaryGroupFromGeometry();

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step / 10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  // Boundary condition (Neumann)
  Matrix<Real> stress = Matrix<Real, 2, 2>::Identity() * 1e3;

  model.applyBC(BC::Neumann::FromHigherDim(stress), "boundary_0");

  model.setBaseName("square-orthotrope");
  model.addDumpFieldVector("displacement");
  model.addDumpField("mass");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.dump();

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for (Int s = 0; s < max_steps; ++s) {
    model.solveStep();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
           << std::endl;

    if (s % 100 == 0) {
      model.dump();
    }
  }

  energy.close();

  finalize();

  return EXIT_SUCCESS;
}
