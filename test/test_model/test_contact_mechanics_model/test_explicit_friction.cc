/**
 * @file   test_explicit_friction.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Sun Jun 06 2021
 * @date last modification: Sun Jun 06 2021
 *
 * @brief  Test contact mechanics with friction
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
template <typename T> std::vector<T> arrange(T start, T stop, T step = 1) {
  std::vector<T> values;
  for (T value = start; value <= stop; value += step)
    values.push_back(value);
  return values;
}

int main(int argc, char * argv[]) {

  UInt max_normal_steps = 2500;
  UInt max_shear_steps = 7500;
  Real max_shear_displacement = 1e-1;
  Real max_normal_displacement = 2e-2;
  Real damping_ratio = 0.99;

  std::string mesh_file = "sliding-block-2D.msh";
  std::string material_file = "material-friction.dat";

  const UInt spatial_dimension = 2;

  initialize(material_file, argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read(mesh_file);

  CouplerSolidContact coupler(mesh);

  auto & solid = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && material_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              solid);
  solid.setMaterialSelector(material_selector);

  coupler.initFull(_analysis_method = _explicit_lumped_mass);

  auto && surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "lower");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "lower");

  Real time_step = solid.getStableTimeStep();
  time_step *= 0.05;

  coupler.setTimeStep(time_step);

  std::cout << "Stable time increment    : " << time_step << " sec "
            << std::endl;

  coupler.setBaseName("explicit-friction");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("tangential_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("strain");
  coupler.addDumpField("stress");
  coupler.addDumpField("contact_state");

  auto & velocity = solid.getVelocity();
  auto & gaps = contact.getGaps();

  auto xi = arrange<Real>(0, 1, 1. / max_shear_steps);

  std::vector<Real> shear_displacements;
  std::transform(xi.begin(), xi.end(), std::back_inserter(shear_displacements),
                 [&](Real & p) -> Real {
                   return 0. + (max_shear_displacement)*pow(p, 3) *
                                   (10 - 15 * p + 6 * pow(p, 2));
                 });

  auto normal_xi = arrange<Real>(0, 1, 1. / max_normal_steps);

  std::vector<Real> normal_displacements;
  std::transform(normal_xi.begin(), normal_xi.end(),
                 std::back_inserter(normal_displacements),
                 [&](Real & p) -> Real {
                   return 0. + (max_normal_displacement)*pow(p, 3) *
                                   (10 - 15 * p + 6 * pow(p, 2));
                 });

  auto max_steps = max_normal_steps + max_shear_steps;

  auto & contact_nodes = surface_selector->getSlaveList();
  auto & tangential_traction = contact.getTangentialTractions();

  for (UInt s : arange(max_steps)) {

    if (s < max_normal_steps) {
      solid.applyBC(BC::Dirichlet::FixedValue(-normal_displacements[s], _y),
                    "loading");
    } else {
      solid.applyBC(BC::Dirichlet::FixedValue(
                        shear_displacements[s - max_normal_steps], _x),
                    "loading");
    }

    coupler.solveStep();

    for (auto && tuple : zip(gaps, make_view(velocity, spatial_dimension))) {
      auto & gap = std::get<0>(tuple);
      auto & vel = std::get<1>(tuple);
      if (gap > 0) {
        vel *= damping_ratio;
      }
    }

    if (s % 100 == 0) {
      coupler.dump();
    }

    auto sum = std::accumulate(tangential_traction.begin(),
                               tangential_traction.end(), 0.0);
    auto num_tang_traction = std::abs(sum) / contact_nodes.size();

    Real exp_tang_traction = 0.3 * 1.4e6;

    Real error =
        std::abs(num_tang_traction - exp_tang_traction) / exp_tang_traction;

    if (error > 1e-3 and num_tang_traction > exp_tang_traction) {
      std::cerr << error << "----" << num_tang_traction << std::endl;
      return EXIT_FAILURE;
    }
  }

  coupler.dump();

  finalize();
  return EXIT_SUCCESS;
}
