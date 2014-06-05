/**
 * @file   test_solid_mechanics_model_cube3d_tetra10.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 *
 * @date   Mon Dec 06 14:43:53 2010
 *
 * @brief  test of the class SolidMechanicsModel on the 3d cube
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize("material.dat", argc, argv);
  akantu::UInt max_steps = 1000;
  akantu::Real epot, ekin;

  akantu::Mesh mesh(3);

  mesh.read("cube2.msh");

  akantu::SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  akantu::Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;


  /// boundary conditions
  akantu::Real eps = 1e-2;
  akantu::UInt nb_nodes = model.getFEEngine().getMesh().getNbNodes();
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    model.getDisplacement().storage()[3*i] = model.getFEEngine().getMesh().getNodes().storage()[3*i] / 100.;

    if(model.getFEEngine().getMesh().getNodes().storage()[3*i] <= eps) {
      model.getBlockedDOFs().storage()[3*i    ] = true;
    }

    if(model.getFEEngine().getMesh().getNodes().storage()[3*i + 1] <= eps) {
      model.getBlockedDOFs().storage()[3*i + 1] = true;
    }
  }
  //  model.getDisplacement().storage()[1] = 0.1;
  model.setBaseName("cube3d_t10");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  std::ofstream energy;
  energy.open("energy.csv");
  energy << "id,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();


    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    std::cerr << "passing step " << s << "/" << max_steps << std::endl;
    energy << s << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

    if(s % 10 == 0) model.dump();
  }

  energy.close();
  akantu::finalize();
  return EXIT_SUCCESS;
}
