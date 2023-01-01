/**
 * @file   contact_mechanics_internodes_model.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Thu Jul 17 2022
 *
 * @brief Model for Contact Mechanics Internodes
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
#include "contact_mechanics_internodes_model.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
ContactMechanicsInternodesModel::ContactMechanicsInternodesModel(
    Mesh & mesh, UInt dim, const ID & id,
    std::shared_ptr<DOFManager> dof_manager, const ModelType model_type)
    : Model(mesh, model_type, dof_manager, dim, id) {

  this->detector = std::make_unique<ContactDetectorInternodes>(
      this->mesh, id + ":contact_detector");

  this->solid = std::make_unique<SolidMechanicsModel>(
      this->mesh, spatial_dimension, id + ":solid_mechanics_model",
      this->dof_manager);
}

/* -------------------------------------------------------------------------- */
ContactMechanicsInternodesModel::~ContactMechanicsInternodesModel() = default;

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::initFullImpl(
    const ModelOptions & options) {
  // check if static solver
  if (options.analysis_method != _static) {
    AKANTU_EXCEPTION(options.analysis_method
        << " is not a valid analysis method for "
        "ContactMechanicsInternodesModel");
  }

  // run contact detection
  detector->findContactNodes(detector->getMasterNodeGroup(), detector->getSlaveNodeGroup());

  Model::initFullImpl(options);

  // init solid mechanics model
  solid->initFull(options);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleResidual() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.assembleResidual();

  auto & master_node_group = detector->getMasterNodeGroup();
  auto & initial_master_node_group = detector->getInitialMasterNodeGroup();
  auto & slave_node_group = detector->getSlaveNodeGroup();

  Vector<Real> positions_master(master_node_group.size()*spatial_dimension);
  Vector<Real> positions_slave(slave_node_group.size()*spatial_dimension);

  for (auto && node_data : enumerate(master_node_group.getNodes())) {
    auto i = std::get<0>(node_data);
    auto node = std::get<1>(node_data);

    auto pos = detector->getNodePosition(node);
    for (int dim = 0; dim < spatial_dimension; dim++) {
      positions_master(spatial_dimension*i + dim) = pos(dim);
    }
  }

  for (auto && node_data : enumerate(slave_node_group.getNodes())) {
    auto i = std::get<0>(node_data);
    auto node = std::get<1>(node_data);

    auto pos = detector->getNodePosition(node);
    for (int dim = 0; dim < spatial_dimension; dim++) {
      positions_slave(spatial_dimension*i + dim) = pos(dim);
    }
  }

  // interpolation matrix
  auto && R_master_slave_ext = detector->constructInterpolationMatrix(
      master_node_group, slave_node_group, detector->getSlaveRadiuses());

  // calculate differences
  Array<Real> differences(initial_master_node_group.size(),
      spatial_dimension, 0.);

  // scale with E for lower condition number
  // TODO: correct lambdas with lambda / E
  Real E = this->solid->getMaterial(0).get("E");

  // R_master_slave_ext * positions_slave - positions_master
  Vector<Real> tmp = (R_master_slave_ext * positions_slave - positions_master) * E;

  UInt next_active_node = 0;
  for (auto && node_data : enumerate(initial_master_node_group.getNodes())) {
    auto i = std::get<0>(node_data);
    auto node = std::get<1>(node_data);

    if (master_node_group.find(node) != -1) {
      UInt active_node_index = next_active_node++;
      for (int dim = 0; dim < spatial_dimension; dim++) {
        differences(i, dim) = tmp(spatial_dimension*active_node_index + dim);
      }
    }
  }

  this->dof_manager->assembleToResidual("lambdas", differences, 1);
}

/* -------------------------------------------------------------------------- */
MatrixType
ContactMechanicsInternodesModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _unsymmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    assembleInternodesMatrix();
  }

  if (matrix_id == "M") {
    solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleLumpedMatrix(
    const ID & /*matrix_id*/) {

  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleInternodesMatrix() {
  solid->assembleStiffnessMatrix();

  auto & initial_master_node_group = detector->getInitialMasterNodeGroup();
  auto & initial_slave_node_group = detector->getInitialSlaveNodeGroup();

  auto & master_node_group = detector->getMasterNodeGroup();
  auto & slave_node_group = detector->getSlaveNodeGroup();

  // R matrices, need for calculations of blocks
  auto && R_master_slave = detector->constructInterpolationMatrix(
      master_node_group, slave_node_group, detector->getSlaveRadiuses());
  auto && R_slave_master = detector->constructInterpolationMatrix(
      slave_node_group, master_node_group, detector->getMasterRadiuses());

  // get interface masses
  auto && M_master = assembleInterfaceMass(master_node_group);
  auto && M_slave = assembleInterfaceMass(slave_node_group);

  // assemble B and B_tilde
  // B_master = - M_master;
  auto B_slave = M_slave * R_slave_master;

  //B_tilde_master = I;
  //B_tilde_slave = - R_master_slave;

  TermsToAssemble termsB("displacement", "lambdas");
  TermsToAssemble termsB_tilde("lambdas", "displacement");

  // Only used with a "diagonal 1" for "inactive" nodes, to make the matrix non-singular
  TermsToAssemble termsLambdaDiag("lambdas", "lambdas");

  Real E = this->solid->getMaterial(0).get("E");

  // fill in active nodes
  for (auto && ref_node_data :
       enumerate(initial_master_node_group.getNodes())) {
    auto i = std::get<0>(ref_node_data);
    auto ref_node = std::get<1>(ref_node_data);

    auto local_i = master_node_group.find(ref_node);

    if (local_i != -1) {
      for (auto && master_node_data :
           enumerate(master_node_group.getNodes())) {
        auto j = std::get<0>(master_node_data);
        auto master_node = std::get<1>(master_node_data);

        for (int dim = 0; dim < spatial_dimension; dim++) {
          UInt idx_i = local_i*spatial_dimension + dim;
          UInt idx_j = j*spatial_dimension + dim;
          UInt global_idx_i = i*spatial_dimension + dim;
          UInt global_idx_j = master_node*spatial_dimension + dim;

          termsB(global_idx_j, global_idx_i) = - E*M_master(idx_j, idx_i);

          if (idx_i == idx_j) {
            termsB_tilde(global_idx_i, global_idx_j) = E;
          }
        }
      }
    } else {
      for (int s = 0; s < spatial_dimension; s++) {
        UInt global_idx = i * spatial_dimension + s;
        termsLambdaDiag(global_idx, global_idx) = E;
      }
    }
  }

  for (auto && ref_node_data :
       enumerate(initial_master_node_group.getNodes())) {
    auto i = std::get<0>(ref_node_data);
    auto ref_node = std::get<1>(ref_node_data);

    auto local_i = master_node_group.find(ref_node);

    if (local_i != UInt(-1)) {
      for (auto && slave_node_data :
           enumerate(slave_node_group.getNodes())) {
        auto j = std::get<0>(slave_node_data);
        auto slave_node = std::get<1>(slave_node_data);

        for (int dim = 0; dim < spatial_dimension; dim++) {
          UInt idx_i = local_i*spatial_dimension + dim;
          UInt idx_j = j*spatial_dimension + dim;
          UInt global_idx_i = i*spatial_dimension + dim;
          UInt global_idx_j = slave_node*spatial_dimension + dim;

          termsB(global_idx_j, global_idx_i) = E * B_slave(idx_j, idx_i);

          termsB_tilde(global_idx_i, global_idx_j) =
              - E * R_master_slave(idx_i, idx_j);
        }
      }
    } else {
      for (int s = 0; s < spatial_dimension; s++) {
        UInt global_idx = i * spatial_dimension + s;
        termsLambdaDiag(global_idx, global_idx) = E;
      }
    }
  }

  this->dof_manager->assemblePreassembledMatrix("K", termsB);
  this->dof_manager->assemblePreassembledMatrix("K", termsB_tilde);
  this->dof_manager->assemblePreassembledMatrix("K", termsLambdaDiag);
}

/* -------------------------------------------------------------------------- */
Matrix<Real> ContactMechanicsInternodesModel::assembleInterfaceMass(
    const NodeGroup & contact_node_group) {

  // get interface mass for a node group
  this->assembleMatrix("M");
  const auto & M = this->dof_manager->getMatrix("M");

  auto nb_contact_nodes = contact_node_group.size();

  Matrix<Real> M_contact(nb_contact_nodes*spatial_dimension,
      nb_contact_nodes*spatial_dimension);

  for (auto && ref_node_data : enumerate(contact_node_group.getNodes())) {
      auto i = std::get<0>(ref_node_data);
      auto ref_node = std::get<1>(ref_node_data);

    for (auto && node_data :
        enumerate(contact_node_group.getNodes())) {
      auto j = std::get<0>(node_data);
      auto node = std::get<1>(node_data);

      for (int dim = 0; dim < spatial_dimension; dim++) {
        UInt global_idx_ref = ref_node*spatial_dimension + dim;
        UInt global_idx = node*spatial_dimension + dim;

        UInt idx_ref = i*spatial_dimension + dim;
        UInt idx = j*spatial_dimension + dim;

        // M must be constant otherwise the modifying overload of operator() is
        // used, which throws if the entry is absent
        M_contact(idx_ref, idx) = M(global_idx_ref, global_idx);
      }
    }
  }

  return M_contact;
}

std::set<UInt> ContactMechanicsInternodesModel::findPenetratingNodes(
    const NodeGroup & ref_group, const NodeGroup & eval_group, const Array<Real> & eval_radiuses) {
  // TODO don't hardcode this, should be mesh_size * relative_tolerance
  const Real penetration_tolerance = -0.01;

  const auto R_ref_eval = detector->constructInterpolationMatrix(ref_group, eval_group, eval_radiuses);

  std::set<UInt> ref_penetration_nodes;
  for (auto && entry : enumerate(ref_group)) {
    auto i = std::get<0>(entry);
    auto ref_node = std::get<1>(entry);

    // 1) Compute gap: ref gaps = R_ref_eval * eval_positions - ref_positions
    // (pointing towards the inside of the ref body if penetrating)
    Vector<Real> gap(spatial_dimension);
    for (auto && inner_entry : enumerate(eval_group)) {
      auto j = std::get<0>(inner_entry);
      auto eval_node = std::get<1>(inner_entry);

      for (UInt s : arange(spatial_dimension)) {
        gap(s) += R_ref_eval(i * spatial_dimension + s, j * spatial_dimension + s)
                  * detector->getPositions()(eval_node, s);
      }
    }
    gap -= detector->getNodePosition(ref_node);

    // 2) Compute normal at ref_node by averaging element normals
    auto normal = getInterfaceNormalAtNode(ref_group, ref_node);

    // 3) Penetration if dot(gap, normal) is negative (with some tolerance)
    if (normal.dot(gap) < penetration_tolerance) {
      ref_penetration_nodes.insert(ref_node);
    }
  }

  return ref_penetration_nodes;
}

/* -------------------------------------------------------------------------- */
bool ContactMechanicsInternodesModel::updateAfterStep() {
  // Update positions
  detector->getPositions().copy(solid->getCurrentPosition());

  auto & initial_master_node_group = detector->getInitialMasterNodeGroup();
  auto & master_node_group = detector->getMasterNodeGroup();
  auto & slave_node_group = detector->getSlaveNodeGroup();

  // Find nodes to remove (if their projected Lagrange multiplier is positive)
  std::set<UInt> to_remove_master, to_remove_slave;
  {
    // Master
    // The lambdas are actually lambdas/E, but that's fine for our sign check
    const auto & lambda_solution = dof_manager->getSolution("lambdas");
    Vector<Real> lambdas_master(master_node_group.size() * spatial_dimension);
    // We need to filter out unused lambdas
    UInt next_used_lambda = 0;
    for (auto && entry : enumerate(initial_master_node_group)) {
      UInt i = std::get<0>(entry);
      UInt maybe_active_node = std::get<1>(entry);
      if (master_node_group.find(maybe_active_node) != -1) {
        for (UInt s = 0; s < spatial_dimension; ++s) {
          lambdas_master(next_used_lambda++) =
              lambda_solution(i * spatial_dimension + s);
        }
      }
    }

    for (auto && entry : enumerate(master_node_group)) {
      auto i = std::get<0>(entry);
      auto node = std::get<1>(entry);

      auto normal = getInterfaceNormalAtNode(master_node_group, node);
      Real dot_product = 0;
      for (UInt s : arange(spatial_dimension)) {
        dot_product += lambdas_master(i * spatial_dimension + s) * normal(s);
      }

      if (dot_product > 0) {
        to_remove_master.insert(node);
      }
    }

    // Slave
    // Interpolate the Lagrange multipliers of the slave
    auto old_R21 = detector->constructInterpolationMatrix(
        slave_node_group, master_node_group, detector->getMasterRadiuses());

    auto lambdas_slave = old_R21 * lambdas_master;
    lambdas_slave *= -1;

    for (auto && entry : enumerate(slave_node_group)) {
      auto i = std::get<0>(entry);
      auto node = std::get<1>(entry);

      auto normal = getInterfaceNormalAtNode(slave_node_group, node);
      Real dot_product = 0;
      for (UInt s : arange(spatial_dimension)) {
        dot_product += lambdas_slave(i * spatial_dimension + s) * normal(s);
      }

      if (dot_product > 0) {
        to_remove_slave.insert(node);
      }
    }
  }

  // Find nodes to add (if they are penetrating)
  std::set<UInt> to_add_master, to_add_slave;
  {
    // COMPUTE PENETRATION
    // Start with all the initial nodes
    NodeGroup penetration_master_group("penetration_master", mesh);
    NodeGroup penetration_slave_group("penetration_slave", mesh);
    penetration_master_group.append(detector->getInitialMasterNodeGroup());
    penetration_slave_group.append(detector->getInitialSlaveNodeGroup());

    // Find contact nodes and radii with the contact detection algorithm.
    // This "pollutes" the state of the detector, but it doesn't matter because we won't be using it again in this iteration.
    detector->findContactNodes(penetration_master_group, penetration_slave_group);

    // Compute penetrating nodes (i.e. nodes to add to the interface)
    to_add_master =
        findPenetratingNodes(penetration_master_group,
                                       penetration_slave_group,
                                       detector->getSlaveRadiuses());
    to_add_slave =
        findPenetratingNodes(penetration_slave_group,
                                       penetration_master_group,
                                       detector->getMasterRadiuses());
  }

  // Remove nodes with positive projected Lagrange multipliers
  bool interface_nodes_changed = false;
  if (master_node_group.applyNodeFilter([&] (UInt master_node) {
        bool keep = to_remove_master.count(master_node) == 0;
        return keep;
      }) > 0) {
    interface_nodes_changed = true;
  }
  if (slave_node_group.applyNodeFilter([&] (UInt slave_node) {
        bool keep = to_remove_slave.count(slave_node) == 0;
        return keep;
      }) > 0) {
    interface_nodes_changed = true;
  }

  // Add nodes that are penetrating
  for (auto node : master_node_group) {
    to_add_master.erase(node);
  }
  for (auto node : slave_node_group) {
    to_add_slave.erase(node);
  }

  for (auto node : to_add_master) {
    master_node_group.add(node);
    interface_nodes_changed = true;
  }
  for (auto node : to_add_slave) {
    slave_node_group.add(node);
    interface_nodes_changed = true;
  }

  if (interface_nodes_changed) {
    // Make sure to re-optimize the node groups after a modification!
    master_node_group.optimize();
    slave_node_group.optimize();
  }

  return interface_nodes_changed;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::solveStep(SolverCallback & callback,
    const ID & solver_id) {
  ModelSolver::solveStep(callback, solver_id);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::solveStep(const ID & solver_id) {

  ModelSolver::solveStep(solver_id);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::predictor() {
  auto & solid_callback = aka::as_type<SolverCallback>(*solid);
  solid_callback.predictor();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::corrector() {
  auto & solid_callback = aka::as_type<SolverCallback>(*solid);
  solid_callback.corrector();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::beforeSolveStep() {
  auto & solid_callback = aka::as_type<SolverCallback>(*solid);
  solid_callback.beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::afterSolveStep(bool converged) {
  auto & solid_callback = aka::as_type<SolverCallback>(*solid);
  solid_callback.afterSolveStep(converged);
}

void ContactMechanicsInternodesModel::printself(std::ostream & stream,
                                                int indent) const {
  AKANTU_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType
ContactMechanicsInternodesModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
ContactMechanicsInternodesModel::getDefaultSolverID(
    const AnalysisMethod & method) {
  switch (method) {
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions ContactMechanicsInternodesModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_linear;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type for "
                             "ContactMechanicsInternodesModel");
    break;
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::initSolver(
    TimeStepSolverType time_step_solver_type,
    NonLinearSolverType non_linear_solver_type) {

  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.initSolver(time_step_solver_type, non_linear_solver_type);

  // Note that we always allocate the lambdas for all the initial master nodes.
  // If we decide to remove some initial master nodes from the interface, we
  // will fill the system with an "identity matrix" for the removed nodes.
  auto nb_initial_master_nodes = detector->getInitialMasterNodeGroup().size();

  // allocate lambdas
  this->lambdas = std::make_unique<Array<Real>>(nb_initial_master_nodes,
      spatial_dimension, id + ":lambdas");

  // allocate blocked dofs for lambdas
  this->blocked_dofs = std::make_unique<Array<bool>>(
      nb_initial_master_nodes, spatial_dimension, false, id + ":blocked_dofs");

  this->dof_manager->registerDOFs("lambdas", *lambdas, _dst_generic);
  this->dof_manager->registerBlockedDOFs("lambdas", *blocked_dofs);
}

} // namespace akantu
