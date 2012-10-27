/**
 * @file   test_solid_mechanics_model_cube3d.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Tue Aug 17 11:31:22 2010
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
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::UInt max_steps = 10000;
  akantu::Real epot, ekin;

  akantu::Mesh mesh(3);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("cube_structured.msh", mesh);

  akantu::SolidMechanicsModel model(mesh);

  model.setPBC(1,1,1);
  model.initFull("material.dat");
  akantu::Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// boundary conditions
  akantu::UInt nb_nodes = model.getFEM().getMesh().getNbNodes();
  akantu::Real eps = 1e-16;
  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    model.getDisplacement().values[3*i] = model.getFEM().getMesh().getNodes().values[3*i] / 100.;

    if(model.getFEM().getMesh().getNodes().values[3*i] <= eps) {
      model.getBoundary().values[3*i    ] = true;
    }

    if(model.getFEM().getMesh().getNodes().values[3*i + 1] <= eps) {
      model.getBoundary().values[3*i + 1] = true;
    }

  }

  model.setBaseName("cube3d_pbc");
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
