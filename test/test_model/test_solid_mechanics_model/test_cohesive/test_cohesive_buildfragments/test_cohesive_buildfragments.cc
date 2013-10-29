/**
 * @file   test_cohesive_buildfragments.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date   Wed Jun 13 11:29:49 2012
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
#include "mesh_utils.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"

//#include "io_helper.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize(argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  const UInt max_steps = 200;
  ElementType type_facet = _segment_2;

  Mesh mesh(spatial_dimension);
  mesh.read("mesh.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull("material.dat", _explicit_lumped_mass, true);

  Real time_step = model.getStableTimeStep()*0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  model.assembleMassLumped();

  Real strain_rate = 1.e5;
  Real L = 0.03;
  Real theoretical_mass = L * L/20. * 2500;
  Real disp_increment = strain_rate * L / 2. * time_step;

  Real epsilon = std::numeric_limits<Real>::epsilon();

  Array<Real> & velocity = model.getVelocity();
  Array<bool> & boundary = model.getBoundary();
  Array<Real> & displacement = model.getDisplacement();
  //  const Array<Real> & residual = model.getResidual();

  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  /// initial conditions
  for (UInt n = 0; n < nb_nodes; ++n)
    velocity(n, 0) = strain_rate * position(n, 0);

  /// boundary conditions
  for (UInt n = 0; n < nb_nodes; ++n) {
    if (std::abs((position(n, 0) - L/2) / (L/2)) < epsilon) {
      boundary(n, 0) = true;
      displacement(n, 0) += disp_increment;
    }
    if (std::abs((position(n, 0) - (-1.) * L/2) / ( (-1) * L/2)) < epsilon) {
      boundary(n, 0) = true;
      displacement(n, 0) -= disp_increment;
    }
  }

  model.updateResidual();

  // model.setBaseName("extrinsic");
  // model.addDumpFieldArray("displacement");
  // model.addDumpField("velocity"    );
  // model.addDumpField("acceleration");
  // model.addDumpField("residual"    );
  // model.addDumpField("stress");
  // model.addDumpField("strain");
  // model.dump();

  ElementType type_cohesive = FEM::getCohesiveElementType(type_facet);
  //  UInt cohesive_index = model.getCohesiveIndex();
  UInt cohesive_index = 1;

  UInt nb_quad_per_facet = model.getFEM("FacetsFEM").getNbQuadraturePoints(type_facet);
  MaterialCohesive & mat_cohesive
    = dynamic_cast<MaterialCohesive&>(model.getMaterial(cohesive_index));
  const Array<Real> & damage = mat_cohesive.getDamage(type_cohesive);

  const Array<Real> & fragment_mass = model.getFragmentsMass();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();

    model.checkCohesiveStress();

    model.updateAcceleration();
    model.explicitCorr();

    /// apply boundary conditions
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs((position(n, 0) - L/2) / (L/2)) < epsilon) {
	displacement(n, 0) += disp_increment;
      }
      if (std::abs((position(n, 0) - (-1.) * L/2) / ( (-1) * L/2)) < epsilon) {
	displacement(n, 0) -= disp_increment;
      }
    }

    if(s % 1 == 0) {
      //      model.dump();

      std::cout << "passing step " << s << "/" << max_steps << std::endl;

      model.computeFragmentsData();

      /// check number of fragments
      UInt nb_fragment_num = model.getNbFragment();

      UInt nb_cohesive_elements = mesh.getNbElement(type_cohesive);

      UInt nb_fragment = 1;
      for (UInt el = 0; el < nb_cohesive_elements; ++el) {
	UInt q = 0;
	while (q < nb_quad_per_facet &&
	       std::abs(damage(el * nb_quad_per_facet + q) - 1) < epsilon) ++q;

	if (q == nb_quad_per_facet) {
	  ++nb_fragment;
	}
      }

      if (nb_fragment != nb_fragment_num) {
	std::cout << "The number of fragments is wrong!" << std::endl;
	return EXIT_FAILURE;
      }

      /// check mass computation
      Real total_mass = 0.;
      for (UInt frag = 0; frag < nb_fragment_num; ++frag) {
	total_mass += fragment_mass(frag);
      }

      if (std::abs(theoretical_mass - total_mass) > epsilon * theoretical_mass * 100) {
	std::cout << "The fragments' mass is wrong!" << std::endl;
	return EXIT_FAILURE;
      }

    }
  }

  /// check velocities
  UInt nb_fragment = model.getNbFragment();
  const Array<Real> & fragment_velocity = model.getFragmentsVelocity();

  for (UInt frag = 0; frag < nb_fragment; ++frag) {
    Real vel = ((frag + 0.5) / nb_fragment * 2 - 1) * disp_increment / time_step;

    if (std::abs(fragment_velocity(frag) - vel) > 100) {
      std::cout << "The fragments' velocity is wrong!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();

  std::cout << "OK: test_cohesive_buildfragments was passed!" << std::endl;
  return EXIT_SUCCESS;
}
