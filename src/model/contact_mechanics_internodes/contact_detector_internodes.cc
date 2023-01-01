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
#include "contact_detector_shared.hh"
#include "element_group.hh"
#include "node_group.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactDetectorInternodes::ContactDetectorInternodes(Mesh & mesh, const ID & id)
    : Parsable(ParserType::_contact_detector, id),
      AbstractContactDetector(mesh, mesh.getNodes()) {

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

  // Remember to fill nodesToElements as they are needed for computeElementSizes
  auto surface_dimension = spatial_dimension - 1;
  this->mesh.fillNodesToElements(surface_dimension);

  auto element_sizes = computeElementSizes(arange(mesh.getNbNodes()));
  this->max_element_size = element_sizes.max_size;
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::parseSection(const ParserSection & section) {
  this->id_master_nodes = section.getParameterValue<std::string>("master");
  this->id_slave_nodes = section.getParameterValue<std::string>("slave");

  this->relative_grid_spacing_factor = section.getParameter("grid_spacing", 3.0);
  this->relative_penetration_tolerance = section.getParameter("penetration_tolerance", 0.1);
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
void ContactDetectorInternodes::findContactNodes(NodeGroup & master_node_group, NodeGroup & slave_node_group) {
  // We need to find master->slave and slave->master interpolation nodes and
  // radii that satisfy the conditions from [1], page 51, equations (2) and (3).

  Real grid_spacing = this->max_element_size * this->relative_grid_spacing_factor;

  bool still_isolated_nodes = true;
  while (still_isolated_nodes) {
    auto master_grid = constructGrid(master_node_group, grid_spacing);
    auto slave_grid = constructGrid(slave_node_group, grid_spacing);

    // computeRadiuses will produce a result satisfying equation (2).
    auto && nb_slave_nodes_inside_radius =
        computeRadiuses(master_radiuses, master_node_group, master_grid,
                        slave_node_group, slave_grid);
    auto && nb_master_nodes_inside_radius =
        computeRadiuses(slave_radiuses, slave_node_group, slave_grid,
                        master_node_group, master_grid);

    // Check equation (3): if a node is still isolated, remove it and iterate.
    still_isolated_nodes = false;

    if (master_node_group.applyNodeFilter([&](auto && node) {
      return nb_master_nodes_inside_radius[node] > 0;
    }) > 0) {
      still_isolated_nodes = true;
    }

    if (slave_node_group.applyNodeFilter([&](auto && node) {
      return nb_slave_nodes_inside_radius[node] > 0;
    }) > 0) {
      still_isolated_nodes = true;
    }

    // To keep removals in O(1), removing from a NodeGroup might reorder nodes
    master_node_group.optimize();
    slave_node_group.optimize();
  }

  // Check that equation (3) is satisfied.
  if (DEBUG_VERIFY_INTERPOLATION_CONDITIONS) {
    for (UInt eval_node : slave_node_group) {
      bool ok = false;
      for (auto entry : enumerate(master_node_group)) {
        auto j = std::get<0>(entry);
        auto ref_node = std::get<1>(entry);

        Real distance = computeDistance(ref_node, eval_node);

        if (distance <= 0.95 * master_radiuses(j) + 1e-9) {
          ok = true;
        }
      }
      if (!ok) {
        AKANTU_EXCEPTION("Node " << eval_node << " is not in contact with any master node");
      }
    }

    for (UInt eval_node : master_node_group) {
      bool ok = false;
      for (auto entry : enumerate(slave_node_group)) {
        auto j = std::get<0>(entry);
        auto ref_node = std::get<1>(entry);

        Real distance = computeDistance(ref_node, eval_node);

        if (distance <= 0.95 * slave_radiuses(j) + 1e-9) {
          ok = true;
        }
      }
      if (!ok) {
        AKANTU_EXCEPTION("Node " << eval_node << " is not in contact with any slave node");
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
SpatialGrid<UInt> ContactDetectorInternodes::constructGrid(const akantu::NodeGroup & node_group, akantu::Real spacing) const {
  SpatialGrid<UInt> grid(spatial_dimension);
  auto spacing_vector = Vector<Real>(spatial_dimension, spacing);
  grid.setSpacing(spacing_vector);

  for (UInt node : node_group) {
    grid.insert(node, getNodePosition(node));
  }

  return grid;
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
std::map<UInt, UInt> ContactDetectorInternodes::computeRadiuses(
    Array<Real> & attack_radiuses, const NodeGroup & ref_node_group,
    const SpatialGrid<UInt> & ref_grid, const NodeGroup & eval_node_group,
    const SpatialGrid<UInt> & eval_grid) {
  Real c = 0.5;
  Real C = 0.95;

  std::vector<UInt> temp_nodes;

  attack_radiuses.resize(ref_node_group.size());

  std::map<UInt, UInt> nb_neighboor_nodes_inside_radiuses;
  std::map<UInt, UInt> nb_opposite_nodes_inside_radiuses;

  UInt nb_iter = 0;

  while (nb_iter < MAX_RADIUS_ITERATIONS) {
    // maximum number of support nodes
    UInt f = std::floor(1 / (std::pow(1 - c, 4) * (1 + 4 * c)));

    for (auto && ref_node_data : enumerate(ref_node_group.getNodes())) {
      auto j = std::get<0>(ref_node_data);
      auto ref_node = std::get<1>(ref_node_data);

      // Compute radius of attack, i.e. distance to closest neighbor node with
      // a scaling factor of c.
      // This guarantees that [1], page 51, equation (2) is satisfied.
      Real attack_radius = std::numeric_limits<double>::max();
      for (auto neighboor_node : ref_grid.setToNeighboring(getNodePosition(ref_node), temp_nodes)) {
        if (neighboor_node != ref_node) {
          Real distance = computeDistance(ref_node, neighboor_node);
          attack_radius = std::min(distance / c, attack_radius);
        }
      }
      attack_radiuses(j) = attack_radius;

      // compute number of neighboor nodes inside attack radius
      for (auto neighbor_node : ref_grid.setToNeighboring(getNodePosition(ref_node), temp_nodes)) {
        Real distance = computeDistance(ref_node, neighbor_node);
        if (distance <= attack_radius) {
          nb_neighboor_nodes_inside_radiuses[neighbor_node]++;
        }
      }

      // compute number of opposite nodes inside C * attack_radius
      for (auto opposite_node : eval_grid.setToNeighboring(getNodePosition(ref_node), temp_nodes)) {
        Real distance = computeDistance(ref_node, opposite_node);
        if (distance < C * attack_radius) {
          nb_opposite_nodes_inside_radiuses[opposite_node]++;
        }
      }
    }

    // maximum number of neighboors inside radius
    // aka maximum number of supports
    UInt max_nb_supports = 0;
    for (const auto & entry : nb_neighboor_nodes_inside_radiuses) {
      max_nb_supports = std::max(max_nb_supports, entry.second);
    }

    // Enforce strict diagonal dominance by rows.
    // See [1], page 51, equation (4).
    if (max_nb_supports <= f) {
      // ok!
      break;
    }

    // correct maximum number of support nodes and then iterate again
    c = 0.5 * (1 + c);

    nb_neighboor_nodes_inside_radiuses.clear();
    nb_opposite_nodes_inside_radiuses.clear();

    nb_iter++;
  }

  if (nb_iter == MAX_RADIUS_ITERATIONS) {
    AKANTU_EXCEPTION("Could not find suitable radii, maximum number of iterations (" << nb_iter << ") was exceeded");
  }

  // Check that equation (2) is satisfied.
  if (DEBUG_VERIFY_INTERPOLATION_CONDITIONS) {
    for (UInt ref_node : ref_node_group) {
      for (auto entry : enumerate(ref_node_group)) {
        auto j = std::get<0>(entry);
        auto other_ref_node = std::get<1>(entry);

        if (ref_node == other_ref_node) {
          continue;
        }

        Real distance = computeDistance(ref_node, other_ref_node);

        if (distance < c * attack_radiuses(j) - 1e-9) {
          AKANTU_EXCEPTION("Radius of attack is too small, distance between nodes "
                           << ref_node << " and " << other_ref_node << " is "
                           << distance << " but radius is " << attack_radiuses(j));
        }
      }
    }
  }

  return nb_opposite_nodes_inside_radiuses;
}

/* -------------------------------------------------------------------------- */
void ContactDetectorInternodes::computeDistancesToRefNode(
    UInt & ref_node, const NodeGroup & eval_node_group,
    Array<Real> & out_array) {
  for (auto && eval_node_data : enumerate(eval_node_group)) {
    auto i = std::get<0>(eval_node_data);
    auto eval_node = std::get<1>(eval_node_data);

    out_array(i) = computeDistance(ref_node, eval_node);
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
