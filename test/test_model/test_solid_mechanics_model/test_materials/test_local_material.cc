/**
 * @file   test_local_material.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 * @date   Fri Nov 26 00:17:56 2010
 *
 * @brief  test of the class SolidMechanicsModel
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
#include "fem.hh"
#include "local_material_damage.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  UInt max_steps = 2000;
  Real epot, ekin;

  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("barre_trou.msh", mesh);
  mesh.getBoundary().createBoundariesFromMeshData("physical_names");

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull("");
  model.readCustomMaterial<LocalMaterialDamage>("material.dat", "local_damage");
  model.initMaterials();

  Real time_step = model.getStableTimeStep();
  model.setTimeStep(time_step/10.);

  model.assembleMassLumped();

  std::cout << model << std::endl;

  /// Dirichlet boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0.0, BC::Functor::_x), "Fixed");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, BC::Functor::_y), "Fixed");

  // Boundary condition (Neumann)
  Matrix<Real> stress(2,2);
  stress.eye(3e6);
  model.applyBC(BC::Neumann::FromHigherDim(stress), "Traction");

  model.setBaseName("local_material");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpField("damage"      );
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );
  model.dump();

  for(UInt s = 0; s < max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();

    if(s % 100 == 0) std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
			       << std::endl;

    if(s % 100 == 0) model.dump();
  }

  akantu::finalize();
  return EXIT_SUCCESS;
}
