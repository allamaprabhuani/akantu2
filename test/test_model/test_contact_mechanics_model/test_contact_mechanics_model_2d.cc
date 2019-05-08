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
#include "contact_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);
  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("contact_2d.msh");

  ContactMechanicsModel model(mesh);
  model.initFull(_analysis_method = _implicit_contact);

  model.setBaseNameToDumper("contact_mechanics", "contact-2d");
  model.addDumpFieldVectorToDumper("contact_mechanics", "contact_force");
  model.addDumpFieldVectorToDumper("contact_mechanics", "external_force");
  model.addDumpFieldVectorToDumper("contact_mechanics", "normals");
  model.addDumpFieldToDumper("contact_mechanics", "gaps");
  model.addDumpFieldToDumper("contact_mechanics", "areas");

  model.search();
  model.dump("contact_mechanics");
  
  return 0;
}

