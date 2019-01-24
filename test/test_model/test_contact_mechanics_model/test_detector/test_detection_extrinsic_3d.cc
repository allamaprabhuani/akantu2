/**
 * @file   test_contact_detection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Dec 18 2018
 * @date last modification: Wed Dec 18 2018
 *
 * @brief  Test for extrinsic detection 3D
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
#include "aka_common.hh"
#include "contact_detector.hh"
#include "contact_element.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real getAnalyticalGap(Mesh &, UInt &);

const Real radius = 0.5;
const UInt spatial_dimension = 3;

using Elements = std::vector<ContactElement>;

int main(int argc, char *argv[]) {

  /*initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("extrinsic_3d.msh");
  
  Elements contact_elements;
  
  ContactDetector detector(mesh);
  detector.search(contact_elements);

  Real epsilon = 1e-5;
  for (auto & element: contact_elements) {
    auto node = element.getSlave();
    auto gap  = element.getGap();
    auto normal = element.getNormal();

    auto analytical_gap = getAnalyticalGap(mesh, node);
    std::cerr << "analytical  = " << analytical_gap << "  computed = " << gap << std::endl;  
    std::cerr << " normal  = " << normal << std::endl;
    //if (abs(analytical_gap - gap) <= epsilon) {
    //  return EXIT_FAILURE;
    //} 
    }*/

  
  return EXIT_SUCCESS;
}


Real getAnalyticalGap(Mesh & mesh, UInt & node) {
  Vector<Real> pos(spatial_dimension);
  Vector<Real> center(spatial_dimension);

  center(0) = 0.0;
  center(1) = 0.501;
  center(2) = 0.0;
  
  
  auto & positions = mesh.getNodes();
  for (UInt s: arange(spatial_dimension)) {
    pos(s) = positions(node, s);
  }

  Real distance = Math::distance_3d(pos.storage(), center.storage());
  Real analytical_gap = distance - radius;

  return analytical_gap;
}
