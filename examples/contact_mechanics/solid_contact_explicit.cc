/**
 * @file   contact_mechanics_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 24 2019
 * @date last modification: Mon Jan 24 2019
 *
 * @brief  contact mechanics model with penalty resolution
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
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "contact_mechanics_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  
  initialize("material.dat", argc, argv);
  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("hertz_2d.msh");

  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method = _static);

  /*solid.setBaseName("static");
  solid.addDumpFieldVector("displacement");
  solid.addDumpField("blocked_dofs");
  solid.addDumpField("external_force");
  solid.addDumpField("internal_force");*/
    
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top");
  solid.applyBC(BC::Dirichlet::IncrementValue(-0.01, _y), "top");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  /*
  auto & solver = solid.getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-12);
  solver.set("convergence_type", _scc_residual);
  */
  
  solid.solveStep();
  //solid.dump();

  auto current_position = solid.getCurrentPosition();

  ContactMechanicsModel contact(mesh, current_position);
  contact.initFull(_analysis_method = _explicit_contact);

  contact.setBaseName("contact-penalty");
  contact.addDumpField("contact_force");
  contact.addDumpField("blocked_dofs");
  contact.addDumpField("gaps");
   
  contact.search();   
  contact.assembleInternalForces();
  
  //contact.dump();

  contact.setBaseNameToDumper("contact_mechanics", "solid");
  contact.addDumpFieldVectorToDumper("contact_mechanics", "displacement");
  contact.addDumpFieldToDumper("contact_mechanics", "external_force");
  
  Array<Real> & contact_force = contact.getInternalForce();
  Array<Real> & external_force = solid.getExternalForce();

  for (auto && values: zip(make_view(external_force),
			   make_view(contact_force)) ) {
    auto & ext_f = std::get<0>(values);
    auto & cont_f = std::get<1>(values);

    ext_f = cont_f;
  }

  for (auto && values: zip(make_view(external_force),
			   make_view(contact_force)) ) {
    auto ext_f = std::get<0>(values);
    auto cont_f = std::get<1>(values);
    std::cerr << ext_f << "  " << cont_f << std::endl;
  }

  //solid.solveStep();
  contact.dump();
  
  finalize();
  return EXIT_SUCCESS;
}

