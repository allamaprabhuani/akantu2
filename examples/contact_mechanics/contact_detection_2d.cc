/**
 * @file   test_contact_detection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Mon Oct 22 2018
 * @date last modification: Mon Oct 22 2018
 *
 * @brief  Test for contact detection
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
#include "contact_detector.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char* argv[]) {

  initialize("material.dat", argc, argv);
  const UInt spatial_diemnsion = 2;

  Mesh mesh(spatial_diemnsion);
  mesh.read("hertz_2d.msh");

  ContactDetector detector(mesh, "bot_body", "top_body");
  detector.setMasterSurface("rigid");
  detector.setSlaveSurface("contact_surface");
  
  detector.search();
  
  finalize();
}
