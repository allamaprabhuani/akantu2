/**
 * @file   test_explicit_dynamic.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Dec 11 2020
 * @date last modification: Sun Jun 06 2021
 *
 * @brief  Test for dynamic explicit contact
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

  UInt max_steps = 2000;
  Real max_displacement = 1e-2;
  Real damping_ratio = 0.99;

  std::string mesh_file = "flat_on_flat.msh";
  std::string material_file = "material.dat";

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

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "upper");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "lower");

  Real time_step = solid.getStableTimeStep();
  time_step *= 0.05;

  coupler.setTimeStep(time_step);

  std::cout << "Stable time increment    : " << time_step << " sec "
            << std::endl;

  coupler.setBaseName("explicit-dynamic");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("strain");
  coupler.addDumpField("stress");

  auto & velocity = solid.getVelocity();
  auto & gaps = contact.getGaps();

  auto xi = arrange<Real>(0, 1, 1. / max_steps);

  std::vector<Real> displacements;
  std::transform(xi.begin(), xi.end(), std::back_inserter(displacements),
                 [&](Real & p) -> Real {
                   return 0. + (max_displacement)*pow(p, 3) *
                                   (10 - 15 * p + 6 * pow(p, 2));
                 });

  for (UInt s : arange(max_steps)) {

    solid.applyBC(BC::Dirichlet::FixedValue(-displacements[s], _y), "loading");
    solid.applyBC(BC::Dirichlet::FixedValue(displacements[s], _y), "fixed");

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
  }

  coupler.dump();

  const ElementType element_type = _quadrangle_4;
  const Array<Real> & stress_vect =
      solid.getMaterial("upper").getStress(element_type);

  auto stress_it = stress_vect.begin(spatial_dimension, spatial_dimension);
  auto stress_end = stress_vect.end(spatial_dimension, spatial_dimension);

  Real stress_tolerance = 1e-2;

  Matrix<Real> presc_stress{{0, 0}, {0, 7e5}};

  for (; stress_it != stress_end; ++stress_it) {
    const auto & stress = *stress_it;

    Real stress_error =
        (std::abs(stress(1, 1)) - presc_stress(1, 1)) / (presc_stress(1, 1));

    // if error is more than 1%
    if (std::abs(stress_error) > stress_tolerance) {
      std::cerr << "stress error: " << stress_error << " > " << stress_tolerance
                << std::endl;
      std::cerr << "stress: " << stress << std::endl
                << "prescribed stress: " << presc_stress << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();
  return EXIT_SUCCESS;
}
