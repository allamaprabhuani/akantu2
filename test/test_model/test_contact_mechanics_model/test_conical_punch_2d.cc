/**
 * @file   test_contact_mechanics_model.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue May 28 2019
 * @date last modification: Tue May 28 2019
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
  initialize("material_punch.dat", argc, argv);

  auto increment = 1e-3;
  auto nb_steps = 10;
  
  Mesh mesh(spatial_dimension);
  mesh.read("conical_punch_2d.msh");
  
  CouplerSolidContact coupler(mesh);

  auto & solid   = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  auto && selector = std::make_shared<MeshDataMaterialSelector<std::string>>(
		     "physical_names",solid);
  solid.setMaterialSelector(selector);
  
  solid.initFull(  _analysis_method = _static);
  contact.initFull(_analysis_method = _explicit_contact);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top_body");

  coupler.initFull(_analysis_method = _explicit_contact);
  
  auto & solver = coupler.getNonLinearSolver();
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", _scc_residual);
  
  coupler.setBaseName("test-conical-punch-2d");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpFieldVector("normals");
  coupler.addDumpFieldVector("tangents");
  coupler.addDumpFieldVector("contact_force");
  coupler.addDumpFieldVector("external_force");
  coupler.addDumpFieldVector("internal_force");
  coupler.addDumpField("gaps");
  coupler.addDumpField("areas");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");
  
  for (auto i : arange(nb_steps)) {

    std::cerr << "Step " << i << std::endl;  
    solid.applyBC(BC::Dirichlet::IncrementValue(-increment, _y), "top_body");

    coupler.solveStep();       
    coupler.dump();
  }

  finalize();
  return EXIT_SUCCESS;
}

