/**
 * @file   contact_mechanics_penalty.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Jan 21 2019
 * @date last modification: Mon Jan 21 2019
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

int main(int argc, char *argv[]){
  
  initialize("material.dat", argc, argv);
  const UInt spatial_dimension = 3;

  Mesh mesh(spatial_dimension);
  mesh.read("hertz_3d.msh");
  
  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method = _static);

  solid.setBaseName("static-3d");
  solid.addDumpFieldVector("displacement");
  solid.addDumpField("blocked_dofs");
  solid.addDumpField("external_force");
  solid.addDumpField("internal_force");
  solid.addDumpField("stress");
  solid.addDumpField("grad_u");
  
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top_body");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _z), "top_body");

  solid.applyBC(BC::Dirichlet::IncrementValue(-0.01, _y), "top_body");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom_surface");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom_surface");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _z), "bottom_surface");

  auto & solver = solid.getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_solution);
  
  solid.solveStep();

  auto current_position = solid.getCurrentPosition();

  ContactMechanicsModel contact(mesh, current_position);
  contact.initFull(_analysis_method = _explicit_contact);
  
  contact.setBaseNameToDumper("contact_mechanics", "contact-3d");
  contact.addDumpFieldVectorToDumper("contact_mechanics", "contact_force");
  contact.addDumpFieldToDumper("contact_mechanics", "gaps");

  contact.search();   
  contact.assembleInternalForces();

  contact.dump("paraview_all");
  contact.dump("contact_mechanics");
  
  Array<Real> & contact_force = contact.getInternalForce();
  Array<Real> & external_force = solid.getExternalForce();
  Array<bool> & blocked_dofs = solid.getBlockedDOFs();
  
  for (auto && values: zip(make_view(external_force),
			   make_view(contact_force),
  			   make_view(blocked_dofs)) ) {
    auto & ext_f = std::get<0>(values);
    auto & cont_f = std::get<1>(values);
    auto & blocked = std::get<2>(values);

    if (cont_f <= 0) {
      ext_f = cont_f;
    }
    
  }

  solid.solveStep();
  contact.dump("paraview_all");
  
  finalize();
  return EXIT_SUCCESS;
}

