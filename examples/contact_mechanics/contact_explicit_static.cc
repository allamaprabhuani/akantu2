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
#include "surface_selector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  
  const UInt spatial_dimension = 2;
  initialize("material.dat", argc, argv);
   
  Mesh mesh(spatial_dimension);
  mesh.read("hertz.msh");
  
  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
		     "physical_names",solid);
  solid.setMaterialSelector(selector);
  
  solid.initFull(  _analysis_method = _static);
  contact.initFull(_analysis_method = _implicit_contact);

  auto && surface_selector = std::make_shared<PhysicalSurfaceSelector>(mesh);
  contact.getContactDetector().setSurfaceSelector(surface_selector);
 
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top");

  coupler.initFull(_analysis_method = _implicit_contact);
  
  coupler.setBaseName("contact-explicit-static");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("velocity");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");

  auto max_steps = 100u;

  for (auto _ [[gnu::unused]] : arange(max_steps)) {

    auto increment = 1e-4;
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "top"); 

    coupler.solveStep();    
    coupler.dump();
  }

  finalize();
  return EXIT_SUCCESS;
}

