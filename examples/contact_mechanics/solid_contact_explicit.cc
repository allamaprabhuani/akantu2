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
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  
  initialize("material_implicit.dat", argc, argv);
  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("hertz_implicit_2d.msh");

  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method = _static);

  solid.setBaseName("static");
  solid.addDumpFieldVector("displacement");
  solid.addDumpField("blocked_dofs");
  solid.addDumpField("external_force");
  solid.addDumpField("internal_force");
    
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top_body");
  solid.applyBC(BC::Dirichlet::IncrementValue(-0.001, _y), "top_body");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");

  auto & solver = solid.getNonLinearSolver();
  solver.set("max_iterations", 1000);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", _scc_residual);

  
  auto & solid_stiffness =
    const_cast<SparseMatrix &>(solid.getDOFManager().getNewMatrix("K", _symmetric));
  
  solid.assembleInternalForces();
  solid.assembleStiffnessMatrix();
  //solid.solveStep();
  //solid.dump();

  auto current_position = solid.getCurrentPosition();

  ContactMechanicsModel contact(mesh, current_position);
  contact.initFull(_analysis_method = _implicit_contact);
  
  contact.setBaseNameToDumper("contact_mechanics", "contact");
  contact.addDumpFieldVectorToDumper("contact_mechanics", "contact_force");
  contact.addDumpFieldVectorToDumper("contact_mechanics", "external_force");
  contact.addDumpFieldToDumper("contact_mechanics", "gaps");
  contact.addDumpFieldToDumper("contact_mechanics", "areas");

  
  contact.search();   
  contact.assembleInternalForces();
  contact.assembleStiffnessMatrix();
  
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

    if (!blocked) {
      ext_f = cont_f;
    }
    
  }

  auto & contact_stiffness =
    const_cast<SparseMatrix &>(contact.getDOFManager().getMatrix("K"));
  
  //auto & solid_stiffness =
  //  const_cast<SparseMatrix &>(solid.getDOFManager().getNewMatrix("K", _symmetric));

  solid_stiffness.add(contact_stiffness);
    
  solid.solveStep();
  contact.dump("paraview_all");
  
  finalize();
  return EXIT_SUCCESS;
}

