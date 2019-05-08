/**
 * @file   test_contact_detection.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Dec 18 2018
 * @date last modification: Wed Dec 18 2018
 *
 * @brief  Test for extrinsic detection 2D
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
#include "aka_grid_dynamic.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;


const Real radius = 0.1;
const UInt spatial_dimension = 3;

auto analyticalCurvedSlave(Mesh & mesh, const UInt & node) {

  auto & positions = mesh.getNodes();

  Real analytical_gap = positions(node, 1);

  Vector<Real> normal(spatial_dimension);
  normal[0] = 0.0;
  normal[1] = 1.0;
  normal[2] = 0.0;

  Vector<Real> tangent(spatial_dimension);
  tangent[0] = -1.0;
  tangent[1] = 0.0;
  tangent[2] = 0.0;
  
  return std::make_tuple(analytical_gap, normal, tangent);
}


auto checkCurvedSlave(int argc, char *argv[]) {

  initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("implicit_3d.msh");

  std::map<UInt, ContactElement> contact_map;
   
  ContactDetector detector(mesh);

  detector.setSurfaceId<Surface::slave>("curved");
  detector.setSurfaceId<Surface::master>("flat");

  SpatialGrid<UInt> master_grid(spatial_dimension);
  SpatialGrid<UInt> slave_grid(spatial_dimension);

  detector.globalSearch(slave_grid, master_grid);
  detector.localSearch(slave_grid, master_grid);
  detector.constructContactMap(contact_map);

  for (auto & entry : contact_map) {
    const auto & slave = entry.first;
    const auto & element = entry.second;
    const auto & gap = element.gap;
    const auto & normal = element.normal;
    const auto & tangent = element.tangents;

    Real analytical_gap;
    Vector<Real> analytical_normal, analytical_tangent;
    
    std::tie(analytical_gap, analytical_normal, analytical_tangent)
      = analyticalCurvedSlave(mesh, slave);
    
    Real tolerance = 1e-8;
 
    auto gap_error = std::abs(gap - analytical_gap);
    if (gap_error > tolerance) {
      std::cerr << "gap error: " << gap_error << " > " << tolerance
                << std::endl;
      std::cerr << "gap: " << gap << std::endl
                << "analytical gap: " << analytical_gap << std::endl;
      return EXIT_FAILURE;
    }
    
    auto normal_error = normal - analytical_normal;
    if (std::abs(normal_error[0]) > tolerance or std::abs(normal_error[1]) > tolerance
	or std::abs(normal_error[2]) > tolerance) {
      std::cerr << "normal error: " << normal_error << " > " << tolerance
                << std::endl;
      std::cerr << "normal: " << normal << std::endl
                << "analytical normal: " << analytical_normal << std::endl;
      return EXIT_FAILURE;
    }

    auto tangent_trans = tangent.transpose();
    auto tang_1 = Vector<Real>(tangent_trans(0));
    auto tang_2 = Vector<Real>(tangent_trans(1));

    auto tangent_error_1 = tang_1 - analytical_tangent;
    auto tangent_error_2 = tang_2 - analytical_tangent;

    if (std::abs(tangent_error_1[0]) > tolerance or std::abs(tangent_error_1[1]) > tolerance
	or std::abs(tangent_error_1[2]) > tolerance) {
      std::cerr << "tangent error: " << tangent_error_1 << " > " << tolerance
                << std::endl;
      std::cerr << "tangent: " << tang_1 << std::endl
                << "analytical tangent: " << analytical_tangent << std::endl;
      return EXIT_FAILURE;
    }
    
  }
}


int main(int argc, char *argv[])
{
  checkCurvedSlave(argc, argv);
  return 0;
}
