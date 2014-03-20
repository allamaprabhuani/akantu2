/**
 * @file   test_cohesive_extrinsic.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Tue May 08 13:01:18 2012
 *
 * @brief  Test for cohesive elements
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

#include <limits>
#include <fstream>
#include <iostream>


/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material.hh"
// #if defined(AKANTU_USE_IOHELPER)
// #  include "io_helper.hh"
// #endif
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 1000;

  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat", SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));
  Real time_step = model.getStableTimeStep()*0.05;
  model.setTimeStep(time_step);
  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();


  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  CohesiveElementInserter & inserter = model.getElementInserter();
  inserter.setLimit('y', -0.30, -0.20);
  model.updateAutomaticInsertion();

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  Array<Real> & position = mesh.getNodes();
  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  UInt nb_nodes = mesh.getNbNodes();

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
      boundary(n, 1) = true;

    if (position(n, 0) > 0.99 || position(n, 0) < -0.99)
      boundary(n, 0) = true;
  }

  /// initial conditions
  Real loading_rate = 0.5;
  Real disp_update = loading_rate * time_step;
  for (UInt n = 0; n < nb_nodes; ++n) {
    velocity(n, 1) = loading_rate * position(n, 1);
  }

  model.updateResidual();

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("residual"    );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();

  // std::ofstream edis("edis.txt");
  // std::ofstream erev("erev.txt");

  //  Array<Real> & residual = model.getResidual();

  //  const Array<Real> & stress = model.getMaterial(0).getStress(type);

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    /// update displacement on extreme nodes
    for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
      if (position(n, 1) > 0.99 || position(n, 1) < -0.99)
    	displacement(n, 1) += disp_update * position(n, 1);
    }

    model.checkCohesiveStress();
    model.solveStep();

    if(s % 1 == 0) {
      model.dump();

      std::cout << "passing step " << s << "/" << max_steps << std::endl;
    }


    // Real Ed = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getDissipatedEnergy();
    // Real Er = dynamic_cast<MaterialCohesive&> (model.getMaterial(1)).getReversibleEnergy();

    // edis << s << " "
    // 	 << Ed << std::endl;

    // erev << s << " "
    // 	 << Er << std::endl;

  }

  // edis.close();
  // erev.close();

  //  mesh.write("mesh_final.msh");

  Real Ed = model.getEnergy("dissipated");

  Real Edt = 200 * std::sqrt(2);

  std::cout << Ed << " " << Edt << std::endl;

  if (Ed < Edt * 0.999 || Ed > Edt * 1.001 || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  // for (UInt n = 0; n < position.getSize(); ++n) {
  //   for (UInt s = 0; s < spatial_dimension; ++s) {
  //     position(n, s) += displacement(n, s);
  //   }
  // }


  finalize();

  std::cout << "OK: test_cohesive_extrinsic was passed!" << std::endl;
  return EXIT_SUCCESS;
}
