/**
 * @file   phase_field_notch.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Oct 02 2018
 * @date last modification: Wed Apr 07 2021
 *
 * @brief  Example of phase field model
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
#include "coupler_solid_phasefield.hh"
#include "group_manager.hh"
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using millisecond = std::chrono::duration<double, std::milli>;

const UInt spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
class PhaseFieldElementFilter : public GroupManager::ClusteringFilter {
public:
  PhaseFieldElementFilter(const PhaseFieldModel & model,
                          const Real max_damage = 1.)
      : model(model), is_unbroken(max_damage) {}

  bool operator()(const Element & el) const override {

    const Array<UInt> & mat_indexes =
        model.getPhaseFieldByElement(el.type, el.ghost_type);
    const Array<UInt> & mat_loc_num =
        model.getPhaseFieldLocalNumbering(el.type, el.ghost_type);

    const auto & mat = model.getPhaseField(mat_indexes(el.element));

    UInt el_index = mat_loc_num(el.element);
    UInt nb_quad_per_element =
        model.getFEEngine("PhaseFieldFEEngine")
            .getNbIntegrationPoints(el.type, el.ghost_type);

    const Array<Real> & damage_array = mat.getDamage(el.type, el.ghost_type);

    AKANTU_DEBUG_ASSERT(nb_quad_per_element * el_index < damage_array.size(),
                        "This quadrature point is out of range");

    const Real * element_damage =
        damage_array.storage() + nb_quad_per_element * el_index;

    UInt unbroken_quads = std::count_if(
        element_damage, element_damage + nb_quad_per_element, is_unbroken);

    return (unbroken_quads > 0);
  }

private:
  struct IsUnbrokenFunctor {
    IsUnbrokenFunctor(const Real & max_damage) : max_damage(max_damage) {}
    bool operator()(const Real & x) const { return x > max_damage; }
    const Real max_damage;
  };

  const PhaseFieldModel & model;
  const IsUnbrokenFunctor is_unbroken;
};

int main(int argc, char * argv[]) {

  initialize("material_notch.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square_notch.msh");

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  model.initFull(_analysis_method = _static);
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);

  phase.initFull(_analysis_method = _static);

  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");

  model.setBaseName("phase_notch");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  UInt nbSteps = 1500;
  Real increment = 1e-5;

  auto start_time = clk::now();

  for (UInt s = 1; s < nbSteps; ++s) {

    if (s >= 500) {
      increment = 1.e-6;
    }

    if (s % 200 == 0) {
      constexpr char wheel[] = "/-\\|";
      auto elapsed = clk::now() - start_time;
      auto time_per_step = elapsed / s;
      std::cout << "\r[" << wheel[(s / 10) % 4] << "] " << std::setw(5) << s
                << "/" << nbSteps << " (" << std::setprecision(2) << std::fixed
                << std::setw(8) << millisecond(time_per_step).count()
                << "ms/step - elapsed: " << std::setw(8)
                << second(elapsed).count() << "s - ETA: " << std::setw(8)
                << second((nbSteps - s) * time_per_step).count() << "s)"
                << std::string(' ', 20) << std::flush;
    }
    model.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");

    coupler.solve();

    auto energy = phase.getEnergy();

    if (s % 100 == 0) {
      model.dump();
    }
  }

  // Real damage_limit = 0.08;
  // auto global_nb_clusters =
  //   mesh.createClusters(spatial_dimension, "crack",
  //   PhaseFieldElementFilter(phase, damage_limit));

  //
  // auto nb_fragment = mesh.getNbElementGroups(spatial_dimension);

  // model.dumpGroup("crack");
  //
  // std::cout << global_nb_clusters << std::endl;
  // std::cout << nb_fragment << std::endl;

  finalize();
  return EXIT_SUCCESS;
}
