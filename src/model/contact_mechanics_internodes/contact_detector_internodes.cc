/**
 * @file   contact_detector_internodes.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Thu Jul 17 2022
 *
 * @brief Algorithm to detetect contact nodes
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "contact_detector_internodes.hh"
#include "element_group.hh"
#include "node_group.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactDetectorInternodes::ContactDetectorInternodes(Mesh & mesh, const ID & id)
    : Parsable(ParserType::_contact_detector, id), mesh(mesh) {

  this->spatial_dimension = mesh.getSpatialDimension();

  const Parser & parser = getStaticParser();
  const ParserSection & section =
      *(parser.getSubSections(ParserType::_contact_detector).first);

  this->parseSection(section);

  auto & initial_master_node_group = mesh.createNodeGroup("initial_contact_master_nodes");
  auto & initial_slave_node_group = mesh.createNodeGroup("initial_contact_slave_nodes");

  auto & master_node_group = mesh.createNodeGroup("contact_master_nodes");
  auto & slave_node_group = mesh.createNodeGroup("contact_slave_nodes");

  initial_master_node_group.append(
      mesh.getElementGroup(id_master_nodes).getNodeGroup());
  initial_slave_node_group.append(mesh.getElementGroup(id_slave_nodes).getNodeGroup());

  master_node_group.append(
      mesh.getElementGroup(id_master_nodes).getNodeGroup());
  slave_node_group.append(mesh.getElementGroup(id_slave_nodes).getNodeGroup());
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::parseSection(const ParserSection & section) {
  this->id_master_nodes = section.getParameterValue<std::string>("master");
  this->id_slave_nodes = section.getParameterValue<std::string>("slave");
}

/* -------------------------------------------------------------------------- */
NodeGroup & ContactDetectorInternodes::getInitialMasterNodeGroup() {
  return mesh.getNodeGroup("initial_contact_master_nodes");
}

/* -------------------------------------------------------------------------- */
NodeGroup & ContactDetectorInternodes::getInitialSlaveNodeGroup() {
  return mesh.getNodeGroup("initial_contact_slave_nodes");
}

/* -------------------------------------------------------------------------- */
NodeGroup & ContactDetectorInternodes::getMasterNodeGroup() {
  return mesh.getNodeGroup("contact_master_nodes");
}

/* -------------------------------------------------------------------------- */
NodeGroup & ContactDetectorInternodes::getSlaveNodeGroup() {
  return mesh.getNodeGroup("contact_slave_nodes");
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::findContactNodes() {
  auto & master_node_group = getMasterNodeGroup();
  auto & slave_node_group = getSlaveNodeGroup();
  auto & master_radiuses = getMasterRadiuses();
  auto & slave_radiuses = getSlaveRadiuses();

  bool still_isolated_nodes = true;
  int iteration = 0;
  while(still_isolated_nodes) {
    auto && nb_slave_nodes_inside_radius  = computeRadiuses(master_radiuses,
            master_node_group, slave_node_group);
    auto && nb_master_nodes_inside_radius  = computeRadiuses(slave_radiuses,
        slave_node_group, master_node_group);

    still_isolated_nodes = false;

    UInt i = 0;
    for (UInt master_node : master_node_group.getNodes()) {
      if (nb_master_nodes_inside_radius(i) == 0) {
        master_node_group.remove(master_node);
        still_isolated_nodes = true;
      }
    ++i;
    }

    i = 0;
    for (UInt slave_node : slave_node_group.getNodes()) {
      if (nb_slave_nodes_inside_radius(i) == 0) {
        slave_node_group.remove(slave_node);
        still_isolated_nodes = true;
      }
    ++i;
    }

    master_node_group.optimize();
    slave_node_group.optimize();
  };
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ContactDetectorInternodes::constructInterpolationMatrix(
    const NodeGroup & ref_node_group, const NodeGroup & eval_node_group,
    Array<Real> eval_radiuses) {
  auto && phi_eval_eval =
      constructPhiMatrix(eval_node_group, eval_node_group, eval_radiuses);
  auto && phi_ref_eval =
      constructPhiMatrix(ref_node_group, eval_node_group, eval_radiuses);

  Matrix<Real> phi_eval_eval_inv(eval_node_group.size(),
      eval_node_group.size(), 1.);

  phi_eval_eval_inv.inverse(phi_eval_eval);
  auto && interpol_ref_eval = phi_ref_eval * phi_eval_eval_inv;

  Vector<Real> ones(eval_node_group.size(), 1);
  ones.set(1.0);
  for (UInt i : arange(ref_node_group.size())) {
    Real row_sum = 0.;
    for (UInt k : arange(eval_node_group.size())) {
      row_sum = row_sum + interpol_ref_eval(i, k);
    }
    for (UInt j : arange(eval_node_group.size())) {
      interpol_ref_eval(i, j) = interpol_ref_eval(i, j) / row_sum;
    }
  }

  // extended to 2D
  Matrix<Real> interpol_ref_eval_ext(spatial_dimension*interpol_ref_eval.rows(),
      spatial_dimension*interpol_ref_eval.cols(), 0.);

  for (UInt i : arange(interpol_ref_eval.rows())) {
    for (UInt j : arange(interpol_ref_eval.cols())) {
      for (int dim = 0; dim < spatial_dimension; dim++) {
        interpol_ref_eval_ext(spatial_dimension*i+dim, spatial_dimension*j+dim) = interpol_ref_eval(i, j);
      }
    }
  }

  return interpol_ref_eval_ext;
}

/* -------------------------------------------------------------------------- */
Matrix<Real>
ContactDetectorInternodes::constructPhiMatrix(const NodeGroup & ref_node_group,
    const NodeGroup & eval_node_group, Array<Real> & eval_radiuses) {
  auto && positions = mesh.getNodes();
  auto && positions_view = make_view(positions, spatial_dimension).begin();

  Array<Real> distances(eval_node_group.size());
  Matrix<Real> phi(ref_node_group.size(), eval_node_group.size());
  phi.set(0.);

  for (auto && ref_node_data : enumerate(ref_node_group.getNodes())) {
    auto i = std::get<0>(ref_node_data);
    auto ref_node = std::get<1>(ref_node_data);

    computeDistancesToRefNode(ref_node, eval_node_group, distances);

    for (UInt j : arange(eval_node_group.size())) {
      if (distances(j) <= eval_radiuses(j)) {
        phi(i, j) = computeRadialBasisInterpolation(distances(j),
            eval_radiuses(j));
      }
    }
  }

  return phi;
}

/* -------------------------------------------------------------------------- */
Array<UInt>
ContactDetectorInternodes::computeRadiuses(Array<Real> & attack_radiuses,
    const NodeGroup & ref_node_group, const NodeGroup & eval_node_group) {
  // UInt n = 1; // fixed to 1 neareast neighboor
  Real c = 0.5;  // conditition (2)
  Real C = 0.95; // condition (3)
  // maximum number of support nodes
  UInt f = std::floor(1 / (std::pow(1 - c, 4) * (1 + 4 * c)));
  Real d = 0.05;      // tolerance, for radius of "attack" estimation
  UInt max_iter = 10; // maximum number of iterations

  Array<Real> distances_neighboors(ref_node_group.size());
  Array<Real> distances_opposites(eval_node_group.size());

  attack_radiuses.resize(ref_node_group.size());

  Array<UInt> nb_neighboor_nodes_inside_radiuses(ref_node_group.size(), 1, 0);
  Array<UInt> nb_opposite_nodes_inside_radiuses(eval_node_group.size(), 1, 0);

  UInt nb_iter = 0;
  UInt max_nb_supports = std::numeric_limits<int>::max();

  while (max_nb_supports > f && nb_iter < max_iter) {
    for (auto && ref_node_data : enumerate(ref_node_group.getNodes())) {
      auto j = std::get<0>(ref_node_data);
      auto ref_node = std::get<1>(ref_node_data);

      // distances to all nodes on same surface (aka neighboors)
      computeDistancesToRefNode(ref_node, ref_node_group, distances_neighboors);

      // get distances to all nodes on other surface (aka opposites)
      computeDistancesToRefNode(ref_node, eval_node_group, distances_opposites);

      // minimal distance to neighboors i.e radius of "attack"
      // TODO: could be done in computeDistancesToRefNode method
      Real attack_radius = std::numeric_limits<double>::max();

      // compute radius of attack
      for (auto && neighboor_node_data : enumerate(ref_node_group.getNodes())) {
        auto i = std::get<0>(neighboor_node_data);
        auto neighboor_node = std::get<1>(neighboor_node_data);

        if (neighboor_node != ref_node) {
          if (distances_neighboors(i) <= attack_radius) {
            attack_radius = distances_neighboors(i);
          }
        }
      }

      Real correction_radius =
          std::sqrt(d * d + 0.25 * std::pow(attack_radius, 2));
      correction_radius = std::max<Real>(attack_radius, correction_radius);

      if (correction_radius > attack_radius / c) {
        attack_radius = correction_radius / c;
      } else {
        attack_radius = correction_radius;
      }
      attack_radiuses(j) = attack_radius;

      // compute number of neighboor nodes inside radius
      for (UInt i : arange(ref_node_group.size())) {
        if (distances_neighboors(i) < attack_radius) {
          nb_neighboor_nodes_inside_radiuses(i)++;
        }
      }

      // compute number of opposite nodes inside C*radius
      for (UInt i : arange(eval_node_group.size())) {
        if (distances_opposites(i) < C * attack_radius) {
          nb_opposite_nodes_inside_radiuses(i)++;
        }
      }
    }

    // maximum number of neighboors inside radius
    // aka maximum number of supports
    max_nb_supports = 0;
    for (auto nb_neighboors : nb_neighboor_nodes_inside_radiuses) {
      if (nb_neighboors > max_nb_supports) {
        max_nb_supports = nb_neighboors;
      }
    }

    // correct maximum number of support nodes
    if (max_nb_supports > f) {
      c = 0.5 * (1 + c);
      f = floor(1 / (pow(1 - c, 4) * (1 + 4 * c)));

      nb_neighboor_nodes_inside_radiuses.set(0);
      nb_opposite_nodes_inside_radiuses.set(0);

      nb_iter++;
    }
  }

  return nb_opposite_nodes_inside_radiuses;
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::computeDistancesToRefNode(
    UInt & ref_node, const NodeGroup & eval_node_group,
    Array<Real> & out_array) {
  auto && positions = mesh.getNodes();
  auto && positions_view = make_view(positions, spatial_dimension).begin();
  Vector<Real> ref_pos = positions_view[ref_node];

  for (auto && eval_node_data : enumerate(eval_node_group.getNodes())) {
    auto i = std::get<0>(eval_node_data);
    auto eval_node = std::get<1>(eval_node_data);

    auto && pos_eval = positions_view[eval_node];
    out_array(i) = ref_pos.distance(pos_eval);
  }
}

/* -------------------------------------------------------------------------- */
Real ContactDetectorInternodes::computeRadialBasisInterpolation(
    const Real distance, const Real radius) {
  /// rescaled radial basis function: Wendland
  Real ratio = distance / radius;
  Real phi_of_x = std::pow(1 - ratio, 4) * (1 + 4 * ratio);
  return phi_of_x;
}

} // namespace akantu
