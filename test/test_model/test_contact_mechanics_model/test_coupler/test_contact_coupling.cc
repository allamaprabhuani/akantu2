/**
 * @file   test_contact_coupling.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Thu May 16 2019
 *
 * @brief  Test for contact mechanics model class
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "coupler_solid_contact.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {

  const UInt spatial_dimension = 2;
  initialize("material.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("coupling.msh");

  CouplerSolidContact coupler(mesh);

  auto & solid = coupler.getSolidMechanicsModel();
  auto & contact = coupler.getContactMechanicsModel();

  solid.initFull(_analysis_method = _static);
  contact.initFull(_analysis_method = _implicit_contact);

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "bot_body");
  solid.applyBC(BC::Dirichlet::IncrementValue(0.001, _y), "bot_body");

  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "top");
  solid.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "top");

  coupler.initFull(_analysis_method = _implicit_contact);

  coupler.setBaseName("coupling");
  coupler.addDumpFieldVector("displacement");
  coupler.addDumpField("blocked_dofs");
  coupler.addDumpField("external_force");
  coupler.addDumpField("internal_force");
  coupler.addDumpField("grad_u");
  coupler.addDumpField("stress");

  coupler.solveStep();

  contact.dump();

  return 0;
}
