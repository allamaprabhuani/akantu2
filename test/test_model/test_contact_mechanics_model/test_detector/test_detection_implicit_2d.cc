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
#include <set>
#include <tuple>
 
/* -------------------------------------------------------------------------- */

using namespace akantu;


const Real radius = 0.1;
const UInt spatial_dimension = 2;


auto analyticalCurvedSlave(Mesh & mesh, const UInt & node) {

  auto & positions = mesh.getNodes();

  Real analytical_gap = positions(node, 1);

  Vector<Real> normal(spatial_dimension);
  normal[0] = 0.0;
  normal[1] = 1.0;

  Vector<Real> tangent(spatial_dimension);
  tangent[0] = -1.0;
  tangent[1] = 0.0;
  
  return std::make_tuple(analytical_gap, normal, tangent);
}


auto analyticalCurvedMaster(Mesh & mesh, const UInt & node) {

  auto & positions = mesh.getNodes();
   
  Vector<Real> slave_point(spatial_dimension);
  slave_point[0] = positions(node, 0);
  slave_point[1] = positions(node, 1);

  Real slope   = -radius/slave_point[0];  

  Real sign = slave_point[0] < 0 ? -1.0 : 1.0;
  
  Vector<Real> master_point(spatial_dimension);
  master_point[0] = sign* radius/std::sqrt(1 + slope * slope);
  master_point[1] = slope * master_point[0] + radius;

  auto distance = slave_point - master_point;
  Real analytical_gap = Math::norm(spatial_dimension, distance.storage());

  auto normal = distance.normalize();

  Vector<Real> normal_3d(spatial_dimension + 1);
  normal_3d[0] = normal[0];
  normal_3d[1] = normal[1];
  normal_3d[2] = 0.0;

  Vector<Real> outward_3d(spatial_dimension + 1);
  outward_3d[0] = 0.0;
  outward_3d[1] = 0.0;
  outward_3d[2] = 1.0;

  auto tangent_3d = outward_3d.crossProduct(normal_3d);
  
  Vector<Real> tangent(spatial_dimension);
  tangent[0] = tangent_3d[0];
  tangent[1] = tangent_3d[1];
  
  return std::make_tuple(analytical_gap, normal, tangent);
}




auto checkCurvedSlave(int argc, char *argv[]) {

  initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("implicit_2d.msh");

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
    if (std::abs(normal_error[0]) > tolerance or std::abs(normal_error[1]) > tolerance) {
      std::cerr << "normal error: " << normal_error << " > " << tolerance
                << std::endl;
      std::cerr << "normal: " << normal << std::endl
                << "analytical normal: " << analytical_normal << std::endl;
      return EXIT_FAILURE;
    }

    auto tangent_trans = tangent.transpose();
    auto tang = Vector<Real>(tangent_trans(0));

    
    auto tangent_error = tang - analytical_tangent;
    if (std::abs(tangent_error[0]) > tolerance or std::abs(tangent_error[1]) > tolerance) {
      std::cerr << "tangent error: " << tangent_error << " > " << tolerance
                << std::endl;
      std::cerr << "tangent: " << tang << std::endl
                << "analytical tangent: " << analytical_tangent << std::endl;
      return EXIT_FAILURE;
    }
    
  }

  return EXIT_SUCCESS;
}


auto checkCurvedMaster(int argc, char *argv[]) {

  initialize("options.dat", argc, argv);

  Mesh mesh(spatial_dimension);
  mesh.read("implicit_2d.msh");

  std::map<UInt, ContactElement> contact_map;
   
  ContactDetector detector(mesh);

  detector.setSurfaceId<Surface::slave>("flat");
  detector.setSurfaceId<Surface::master>("curved");
  
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
      = analyticalCurvedMaster(mesh, slave);
    
    Real tolerance = 1e-2;
    
    auto gap_error = std::abs(gap - analytical_gap);
    if (gap_error > tolerance) {
      std::cerr << "slave node: " << slave << std::endl;
      std::cerr << "gap error: " << gap_error << " > " << tolerance
                << std::endl;
      std::cerr << "gap: " << gap << std::endl
                << "analytical gap: " << analytical_gap << std::endl;
      return EXIT_FAILURE;
    }
   
    auto normal_error = normal - analytical_normal;
    if (std::abs(normal_error[0]) > tolerance or std::abs(normal_error[1]) > tolerance) {
      std::cerr << "normal error: " << normal_error << " > " << tolerance
                << std::endl;
      std::cerr << "normal: " << normal << std::endl
                << "analytical normal: " << analytical_normal << std::endl;
      return EXIT_FAILURE;
    }
    
    auto tangent_trans = tangent.transpose();
    auto tang = Vector<Real>(tangent_trans(0));
   
    auto tangent_error = tang - analytical_tangent;
    if (std::abs(tangent_error[0]) > tolerance or std::abs(tangent_error[1]) > tolerance) {
      std::cerr << "tangent error: " << tangent_error << " > " << tolerance
                << std::endl;
      std::cerr << "tangent: " << tang << std::endl
                << "analytical tangent: " << analytical_tangent << std::endl;
      return EXIT_FAILURE;
    }

  }

  return EXIT_SUCCESS;
}


int main(int argc, char *argv[])
{

  //checkCurvedSlave(argc, argv);
  checkCurvedMaster(argc, argv);
  
  return EXIT_SUCCESS;
}




