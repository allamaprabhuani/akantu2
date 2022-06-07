/**
 * @file   test_cohesive_buildfragments.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Sun Oct 19 2014
 * @date last modification:  Thu May 09 2019
 *
 * @brief  Test for cohesive elements
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "fragment_manager.hh"
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

static const bool debug_ = true;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Math::setTolerance(1e-14);

  const Int spatial_dimension = 2;
  const Int max_steps = 200;
  Real strain_rate = 1.e5;
  ElementType type = _quadrangle_4;

  Real L = 0.03;
  Real theoretical_mass = L * L / 20. * 2500;

  ElementType type_facet = Mesh::getFacetType(type);
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  Mesh mesh(spatial_dimension);
  mesh.read("mesh.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(_analysis_method = _explicit_lumped_mass,
                 _is_extrinsic = true);

  Real time_step = model.getStableTimeStep() * 0.05;
  model.setTimeStep(time_step);
  //  std::cout << "Time step: " << time_step << std::endl;

  Real disp_increment = strain_rate * L / 2. * time_step;
  model.assembleMassLumped();

  auto & velocity = model.getVelocity();
  const auto & position = mesh.getNodes();
  auto nb_nodes = mesh.getNbNodes();
  auto nb_elements = mesh.getNbElement(type);

  auto & mesh_facets = mesh.getMeshFacets();

  /// initial conditions
  for (Int n = 0; n < nb_nodes; ++n) {
    velocity(n, 0) = strain_rate * position(n, 0);
  }

  /// boundary conditions
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "Left_side");
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "Right_side");

  auto cohesive_index = 1;

  auto nb_quad_per_facet =
      model.getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);
  auto & mat_cohesive =
      dynamic_cast<MaterialCohesive &>(model.getMaterial(cohesive_index));
  const Array<Real> & damage = mat_cohesive.getDamage(type_cohesive);

  FragmentManager fragment_manager(model, false);

  if (debug_) {
    model.setBaseName("buildfragments");
    model.addDumpFieldVector("displacement");
    model.addDumpField("velocity");
    model.addDumpField("acceleration");
    model.addDumpField("internal_force");
    model.addDumpField("stress");
    model.addDumpField("grad_u");
    model.addDumpFieldToDumper("cohesive elements", "damage");
    model.dump();
  }

  Vector<Int> counts(4);

  /// Main loop
  for (Int s = 1; s <= max_steps; ++s) {
    model.checkCohesiveStress();
    model.solveStep();

    if (debug_) {
      model.dump();
      model.dump("cohesive elements");
    }
    /// apply boundary conditions
    model.applyBC(BC::Dirichlet::IncrementValue(-disp_increment, _x),
                  "Left_side");
    model.applyBC(BC::Dirichlet::IncrementValue(disp_increment, _x),
                  "Right_side");

    const auto & elements_to_facets = mesh_facets.getSubelementToElement();

    for (auto && el : element_range(nb_elements, type)) {
      auto && element_to_facets = elements_to_facets.get(el);
      counts.zero();
      for (auto && facet_data : enumerate(element_to_facets)) {
        const auto & connected_elements =
            mesh_facets.getElementToSubelement(std::get<1>(facet_data));
        counts[std::get<0>(facet_data)] = std::count_if(
            connected_elements.begin(), connected_elements.end(),
            [](auto && element) { return element == ElementNull; });
      }
      std::cout << el << " - " << counts << std::endl;
    }

    if (s % 1 == 0) {
      //      model.dump();
      std::cout << "passing step " << s << "/" << max_steps << std::endl;

      fragment_manager.computeAllData();

      /// check number of fragments
      Int nb_fragment_num = fragment_manager.getNbFragment();

      auto nb_cohesive_elements = mesh.getNbElement(type_cohesive);

      Int nb_fragment = 1;
      for (Int el = 0; el < nb_cohesive_elements; ++el) {
        Int q = 0;
        while (q < nb_quad_per_facet &&
               Math::are_float_equal(damage(el * nb_quad_per_facet + q), 1))
          ++q;

        if (q == nb_quad_per_facet) {
          ++nb_fragment;
        }
      }

      if (nb_fragment != nb_fragment_num) {
        AKANTU_EXCEPTION("The number of fragments is wrong! Got: "
                         << nb_fragment_num << " - expected: " << nb_fragment);
      }

      /// check mass computation
      const auto & fragment_mass = fragment_manager.getMass();
      Vector<Real> zeros = Vector<Real>::Zero(spatial_dimension);
      auto total_mass =
          std::accumulate(fragment_mass.begin(spatial_dimension),
                          fragment_mass.end(spatial_dimension), zeros);

      if (!Math::are_float_equal(theoretical_mass, total_mass(0))) {
        AKANTU_EXCEPTION("The fragments' mass is wrong! Got: "
                         << total_mass(0)
                         << " - expected: " << theoretical_mass);
      }
    }
  }

  model.dump();

  /// check velocities
  Int nb_fragment = fragment_manager.getNbFragment();
  const auto & fragment_velocity = fragment_manager.getVelocity();
  const auto & fragment_center = fragment_manager.getCenterOfMass();

  Real fragment_length = L / nb_fragment;
  Real initial_position = -L / 2. + fragment_length / 2.;

  for (Int frag = 0; frag < nb_fragment; ++frag) {
    Real theoretical_center = initial_position + fragment_length * frag;

    if (!Math::are_float_equal(fragment_center(frag, 0), theoretical_center)) {
      AKANTU_EXCEPTION("The fragments' center is wrong!");
    }

    Real initial_vel = fragment_center(frag, 0) * strain_rate;

    Math::setTolerance(100);

    if (!Math::are_float_equal(fragment_velocity(frag), initial_vel)) {
      AKANTU_EXCEPTION("The fragments' velocity is wrong!");
    }
  }

  std::cout << "OK: test_cohesive_buildfragments was passed!" << std::endl;
  return EXIT_SUCCESS;
}
