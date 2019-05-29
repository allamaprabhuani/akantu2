/**
 * @file   test_contact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Apr 30 2019
 * @date last modification: Tue Apr 30 2019
 *
 * @brief  Test for contact mechanics model class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "contact_mechanics_model.hh"
#include "coupler_solid_contact.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  
  const UInt spatial_dimension = 2;
  initialize("material.dat", argc, argv);

  Real time_step;
  Real time_factor = 0.8;
  UInt max_steps = 5000;
    
  Mesh mesh(spatial_dimension);
  mesh.read("contact_hertz_2d.msh");
  
  SolidMechanicsModel solid(mesh);
  
  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
		     "physical_names", solid);
  solid.setMaterialSelector(selector);
    
  solid.initFull(  _analysis_method = _explicit_lumped_mass);

  time_step = solid.getStableTimeStep();
  std::cout << "Time Step = " << time_step * time_factor << "s (" << time_step
            << "s)" << std::endl;
  solid.setTimeStep(time_step * time_factor);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top");
  
  solid.setBaseName("solid-explicit-dynamic");
  solid.addDumpFieldVector("displacement");
  solid.addDumpFieldVector("velocity");
  solid.addDumpFieldVector("external_force");
  solid.addDumpFieldVector("internal_force");
  solid.addDumpField("blocked_dofs");
  solid.addDumpField("grad_u");
  solid.addDumpField("stress");

  auto & velocity = solid.getVelocity();

  Real damping_ratio = 0.7;
  
  for (auto i : arange(max_steps)) {
    auto increment = time_step * 5e-4;
    
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "top"); 
  
    solid.solveStep();

    Real epot = solid.getEnergy("potential");
    Real ekin = solid.getEnergy("kinetic");

    std::cerr << i << "," << i * increment << "," << epot << "," << ekin << ","
	      << epot + ekin << "," << std::endl;

    for (auto & v : make_view(velocity)) {
      v *= damping_ratio;
    }
    
    solid.dump();
  }

  finalize();
  return EXIT_SUCCESS;
}

