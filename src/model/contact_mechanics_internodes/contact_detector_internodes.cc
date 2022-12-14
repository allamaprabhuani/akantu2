/**
 * @file   contact_detector_internodes.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Wed Dec 14 2022
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

  auto & initial_master_node_group =
      mesh.createNodeGroup("initial_contact_master_nodes");
  auto & initial_slave_node_group =
      mesh.createNodeGroup("initial_contact_slave_nodes");

  auto & master_node_group = mesh.createNodeGroup("contact_master_nodes");
  auto & slave_node_group = mesh.createNodeGroup("contact_slave_nodes");

  initial_master_node_group.append(
      mesh.getElementGroup(id_master_nodes).getNodeGroup());
  initial_slave_node_group.append(
      mesh.getElementGroup(id_slave_nodes).getNodeGroup());

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
void ContactDetectorInternodes::findContactNodes(Real C) {
  auto & master_node_group = getMasterNodeGroup();
  auto & slave_node_group = getSlaveNodeGroup();

  // while there still exist isolated nodes, recompute radius parameters
  bool still_isolated_nodes = true;
  while (still_isolated_nodes) {
    still_isolated_nodes = false;

    computeRadiuses(master_radiuses, master_node_group);
    computeRadiuses(slave_radiuses, slave_node_group);

    for (auto && master_node_data : enumerate(master_node_group.getNodes())) {
      auto master_node_index = std::get<0>(master_node_data);
      auto master_node = std::get<1>(master_node_data);

      Array<Real> distances_to_slave_nodes(slave_node_group.size());
      computeDistancesToRefNode(master_node, slave_node_group,
                                distances_to_slave_nodes);

      UInt counter = 0;
      for (int i = 0; i < slave_node_group.size(); i++)
        if (distances_to_slave_nodes(i) <
            C * master_radiuses(master_node_index))
          counter += 1;

      if (counter == 0) {
        master_node_group.remove(master_node);
        still_isolated_nodes = true;
      }
    }

    for (auto && slave_node_data : enumerate(slave_node_group.getNodes())) {
      auto slave_node_index = std::get<0>(slave_node_data);
      auto slave_node = std::get<1>(slave_node_data);

      Array<Real> distances_to_master_nodes(master_node_group.size());
      computeDistancesToRefNode(slave_node, master_node_group,
                                distances_to_master_nodes);

      UInt counter = 0;
      for (int i = 0; i < master_node_group.size(); i++)
        if (distances_to_master_nodes(i) < C * slave_radiuses(slave_node_index))
          counter += 1;

      if (counter == 0) {
        slave_node_group.remove(slave_node);
        still_isolated_nodes = true;
      }
    }

    master_node_group.optimize();
    slave_node_group.optimize();
  }
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ContactDetectorInternodes::constructInterpolationMatrix(
    const NodeGroup & ref_node_group, const NodeGroup & eval_node_group,
    Array<Real> eval_radiuses) {
  auto && phi_eval_eval =
      constructPhiMatrix(eval_node_group, eval_node_group, eval_radiuses);
  auto && phi_ref_eval =
      constructPhiMatrix(ref_node_group, eval_node_group, eval_radiuses);

  Matrix<Real> phi_eval_eval_inv(eval_node_group.size(), eval_node_group.size(),
                                 1.);

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
  Matrix<Real> interpol_ref_eval_ext(
      spatial_dimension * interpol_ref_eval.rows(),
      spatial_dimension * interpol_ref_eval.cols(), 0.);

  for (UInt i : arange(interpol_ref_eval.rows())) {
    for (UInt j : arange(interpol_ref_eval.cols())) {
      for (int dim = 0; dim < spatial_dimension; dim++) {
        interpol_ref_eval_ext(spatial_dimension * i + dim,
                              spatial_dimension * j + dim) =
            interpol_ref_eval(i, j);
      }
    }
  }

  return interpol_ref_eval_ext;
}

/* -------------------------------------------------------------------------- */
Matrix<Real>
ContactDetectorInternodes::constructPhiMatrix(const NodeGroup & ref_node_group,
                                              const NodeGroup & eval_node_group,
                                              Array<Real> & eval_radiuses) {
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
        phi(i, j) = evaluateRadialBasisFunction(distances(j), eval_radiuses(j));
      }
    }
  }

  return phi;
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::computeRadiuses(
    Array<Real> & radius_parameters, const NodeGroup & ref_node_group, Real c,
    Real C) {

  Array<Real> distances_ref_nodes(ref_node_group.size());
  radius_parameters.resize(ref_node_group.size());

  // array keeping track of how many other nodes which are in a node's support
  Array<UInt> nb_supported_nodes(ref_node_group.size(), 1, 0);
  UInt max_nb_supported_nodes;

  // iteratively increase the parameter c until all conditions are satisfied
  while (true) {
    for (auto && ref_node_data : enumerate(ref_node_group.getNodes())) {
      auto ref_node_index = std::get<0>(ref_node_data);
      auto ref_node = std::get<1>(ref_node_data);

      // compute distance of node to all other nodes
      computeDistancesToRefNode(ref_node, ref_node_group, distances_ref_nodes);

      // compute minimum of distances of node to all other nodes (except itself)
      Real minimum_distance = std::numeric_limits<double>::max();
      for (auto && neighbor_node_data : enumerate(ref_node_group.getNodes())) {
        auto neighbor_node_index = std::get<0>(neighbor_node_data);
        auto neighbor_node = std::get<1>(neighbor_node_data);

        // only account for minimum distance of distinct nodes from ref node
        if (neighbor_node != ref_node)
          if (distances_ref_nodes(neighbor_node_index) <= minimum_distance)
            minimum_distance = distances_ref_nodes(neighbor_node_index);
      }

      // set radius parameter to largest value without violating condition
      radius_parameters(ref_node_index) = minimum_distance / c;

      // update the counters of number of other node's supports it belongs to
      for (UInt i : arange(ref_node_group.size()))
        if (distances_ref_nodes(i) < radius_parameters(ref_node_index))
          nb_supported_nodes(i)++;
    }

    // take maximum of number of other nodes which are in a a node's support
    max_nb_supported_nodes = 0;
    for (auto nb_supports : nb_supported_nodes)
      if (nb_supports > max_nb_supported_nodes)
        max_nb_supported_nodes = nb_supports;

    // stop if condition on number of mutually supported nodes is satisfied
    if (max_nb_supported_nodes < 1 / evaluateRadialBasisFunction(c, 1.0))
      break;

    // increase parameter c if no suitable radius parameters were found
    c = (1.0 + c) / 2.0;
    if (c >= C)
      AKANTU_EXCEPTION("No suitable radius parameters found for given c and C");
    nb_supported_nodes.set(0);
  }
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
Real ContactDetectorInternodes::evaluateRadialBasisFunction(const Real distance,
                                                            const Real radius) {
  /// rescaled radial basis function: Wendland
  Real ratio = distance / radius;
  Real phi_of_x = std::pow(1 - ratio, 4) * (1 + 4 * ratio);
  return phi_of_x;
}

} // namespace akantu
