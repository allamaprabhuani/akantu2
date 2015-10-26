/**
 * @file   test_solid_mechanics_model_square.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Sep 22 2010
 * @date last modification: Fri Apr 04 2014
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  //  akantu::initialize("orthotropic.dat", argc, argv);
  akantu::initialize("orthotropic.dat", argc, argv);
  UInt max_steps = 1000;
  Real epot, ekin;

  Mesh mesh(2);
  mesh.read("square.msh");
  mesh.createBoundaryGroupFromGeometry();

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(Real(1e3));
  model.applyBC(BC::Neumann::FromHigherDim(stress), "boundary_0");

  model.setBaseName       ("square-orthotrope"      );
  model.addDumpFieldVector("displacement");
  model.addDumpField      ("mass"        );
  model.addDumpField      ("velocity"    );
  model.addDumpField      ("acceleration");
  model.addDumpFieldVector("force"       );
  model.addDumpField      ("residual"    );
  model.addDumpField      ("stress"      );
  model.addDumpField      ("grad_u"      );
  model.dump();

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    epot = model.getEnergy("potential");
    ekin = model.getEnergy("kinetic");

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if(s % 100 == 0) model.dump();
  }

  energy.close();

  finalize();

  return EXIT_SUCCESS;
}
