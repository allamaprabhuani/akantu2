/**
 * @file   contact_mechanics_internodes_model.cc
 *
 * @author Moritz Waldleben <moritz.waldleben@epfl.ch>
 *
 * @date creation: Thu Jul 09 2022
 * @date last modification: Thu Jul 09 2022
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


  detector->findContactNodes();
  
  // this doesn't work
  Model::initFullImpl(options);

  // init solid mechanics model
  solid->initFull(options);

  // other idea
  // TODO: change back to initial
  // auto nb_initial_master_nodes = detector->getMasterNodeGroup().getNodes().size();

  // // allocate lambdas
  // this->lambdas = std::make_unique<Array<Real>>(nb_initial_master_nodes, spatial_dimension,
  //                                         id + ":lambdas");

  // // allocate blocked dofs for lambdas
  // this->blocked_dofs = std::make_unique<Array<bool>>(
  //     nb_initial_master_nodes, spatial_dimension, false, id + ":blocked_dofs");

  // this->dof_manager->registerDOFs("lambdas", *lambdas, _dst_generic);
  // this->dof_manager->registerBlockedDOFs("lambdas", *blocked_dofs);

  // Model::initFullImpl(options);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleResidual() {
  // get positions of master and slave nodes
  auto && positions = mesh.getNodes();
  auto && positions_view = make_view(positions, spatial_dimension).begin();

  auto & master_node_group = detector->getMasterNodeGroup();
  auto & initial_master_node_group = detector->getMasterNodeGroup();
  auto & slave_node_group = detector->getSlaveNodeGroup();

  Array<Real> positions_master(master_node_group.size(), spatial_dimension);
  Array<Real> positions_slave(slave_node_group.size(), spatial_dimension);

  UInt i = 0;
  for (auto && node : master_node_group.getNodes()) {
    auto pos = positions_view[node];
    for (int dim = 0; dim < spatial_dimension; dim++) {
      positions_master(i, dim) = pos(dim);
    }
    ++i;
  }

  UInt j = 0;
  for (auto && node : slave_node_group.getNodes()) {
    auto pos = positions_view[node];
    for (int dim = 0; dim < spatial_dimension; dim++) {
      positions_slave(j, dim) = pos(dim);
    }
    ++j;
  }

  // interpolation matrix
  auto && R_master_slave = detector->constructInterpolationMatrix(
      master_node_group, slave_node_group, detector->getSlaveRadiuses());

  // calculate constraints
  Array<Real> constraints(initial_master_node_group.size(), spatial_dimension, 0.);

  for (int dim = 0; dim < spatial_dimension; dim++) {
    auto tmp = R_master_slave * positions_slave(dim);
    for (UInt i = 0; i < master_node_group.size(); i++) {
      constraints(i, dim) = tmp[i] - positions_master(i, dim);
    }
  }
  this->dof_manager->assembleToResidual("lambdas", constraints, 1);
}

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
    solid->assembleStiffnessMatrix();

    auto & initial_master_node_group = detector->getMasterNodeGroup();
    auto & initial_slave_node_group = detector->getSlaveNodeGroup();

    auto nb_initial_master_nodes = initial_master_node_group.getNodes().size();
    auto nb_initial_slave_nodes = initial_slave_node_group.getNodes().size();

    auto & master_node_group = detector->getMasterNodeGroup();
    auto & slave_node_group = detector->getSlaveNodeGroup();

    auto nb_master_nodes = master_node_group.getNodes().size();
    auto nb_slave_nodes = slave_node_group.getNodes().size();

    auto nb_nodes = mesh.getNbNodes();
    auto nb_dofs = spatial_dimension * nb_nodes;
    auto nb_constraint_dofs = spatial_dimension * nb_initial_master_nodes;

    // R matrices, need for calculations of blocks
    auto && R_master_slave = detector->constructInterpolationMatrix(
        master_node_group, slave_node_group, detector->getSlaveRadiuses());
    auto && R_slave_master = detector->constructInterpolationMatrix(
        slave_node_group, master_node_group, detector->getMasterRadiuses());

    // get interface masses
    auto && M_master = constructInterfaceMass(master_node_group);
    auto && M_slave = constructInterfaceMass(slave_node_group);

    // assemble B and B_tilde
    // B_master = - M_master;
    auto B_slave = M_slave * R_slave_master;

    //B_tilde_master = I;
    //B_tilde_slave = - R_master_slave;

    TermsToAssemble termsB("displacement", "lambdas");
    TermsToAssemble termsB_tilde("lambdas", "displacement");

    // fill in active nodes
    UInt i = 0;
    for (auto && node_ref : initial_master_node_group.getNodes()) {
      auto local_i = master_node_group.find(node_ref);
      if (local_i != UInt(-1)) {
        UInt j = 0;
        for (auto && master_node : master_node_group.getNodes()) {
          for (int dim = 0; dim < spatial_dimension; dim++) {
            UInt idx_i = local_i*spatial_dimension + dim;
            UInt idx_j = j*spatial_dimension + dim;
            UInt global_idx_i = i*spatial_dimension + dim;
            UInt global_idx_j = master_node*spatial_dimension + dim;

              termsB(global_idx_j, global_idx_i) = - M_master(idx_j, idx_i);

            if (idx_i == idx_j) {
              termsB_tilde(global_idx_i, global_idx_j) = 1;
            }
          }
          ++j;
        }
      }
      ++i;
    }

    i = 0;
    for (auto && node_ref : initial_master_node_group.getNodes()) {
      auto local_i = master_node_group.find(node_ref);
      if (local_i != UInt(-1)) {
        UInt j = 0;
        for (auto && slave_node : slave_node_group.getNodes()) {
          for (int dim = 0; dim < spatial_dimension; dim++) {
            UInt idx_i = local_i*spatial_dimension + dim;
            UInt idx_j = j*spatial_dimension + dim;
            UInt global_idx_i = i*spatial_dimension + dim;
            UInt global_idx_j = slave_node*spatial_dimension + dim;

              termsB(global_idx_j, global_idx_i) = B_slave(idx_j, idx_i);

              termsB_tilde(global_idx_i, global_idx_j) = - R_master_slave(idx_i, idx_j);
          }
          ++j;
        }
      }
      ++i;
    }

    // assemble C (zeros)
    TermsToAssemble termsC("lambdas", "lambdas");

    this->dof_manager->assemblePreassembledMatrix("K", termsB);
    this->dof_manager->assemblePreassembledMatrix("K", termsB_tilde);
    this->dof_manager->assemblePreassembledMatrix("K", termsC);

    auto & A = this->dof_manager->getMatrix("K");
    A.saveMatrix("A.mtx");
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

void ContactMechanicsInternodesModel::solveStep(SolverCallback & callback,
                                                const ID & solver_id) {
  ModelSolver::solveStep(callback, solver_id);
}

void ContactMechanicsInternodesModel::solveStep(const ID & solver_id) {

  ModelSolver::solveStep(solver_id);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::beforeSolveStep() {
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::afterSolveStep(bool converged) {
  // check for convergence here
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

void ContactMechanicsInternodesModel::initSolver(
    TimeStepSolverType time_step_solver_type,
    NonLinearSolverType non_linear_solver_type) {

  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.initSolver(time_step_solver_type, non_linear_solver_type);

  // TODO: change back to initial
  auto nb_initial_master_nodes = detector->getMasterNodeGroup().getNodes().size();

  // allocate lambdas
  this->lambdas = std::make_unique<Array<Real>>(nb_initial_master_nodes, spatial_dimension,
                                          id + ":lambdas");

  // allocate blocked dofs for lambdas
  this->blocked_dofs = std::make_unique<Array<bool>>(
      nb_initial_master_nodes, spatial_dimension, false, id + ":blocked_dofs");

  this->dof_manager->registerDOFs("lambdas", *lambdas, _dst_generic);
  this->dof_manager->registerBlockedDOFs("lambdas", *blocked_dofs);
}

Matrix<Real> ContactMechanicsInternodesModel::constructInterfaceMass(NodeGroup & contact_node_group) {
  // get interface mass for a node group
  this->assembleMatrix("M");
  auto & M = this->dof_manager->getMatrix("M");

  auto nb_contact_nodes = contact_node_group.size();

  Matrix<Real> M_contact(nb_contact_nodes*spatial_dimension,
      nb_contact_nodes*spatial_dimension);

  UInt i = 0;
  for (auto && node_ref : contact_node_group.getNodes()) {
    UInt j = 0;
    for (auto && node : contact_node_group.getNodes()) {
      for (int dim = 0; dim < spatial_dimension; dim++) {
        UInt global_idx_ref = node_ref*spatial_dimension + dim;
        UInt global_idx = node*spatial_dimension + dim;

        UInt idx_ref = i*spatial_dimension + dim;
        UInt idx = j*spatial_dimension + dim;

        try {
          M_contact(idx_ref, idx) = M(global_idx_ref, global_idx);
        }
        catch(...) {
        M_contact(idx_ref, idx) = 0.;
        }

      }
      ++j;
    }
    ++i;
  }

  return M_contact;
}

void ContactMechanicsInternodesModel::assembleInternodesMatrix() {
  this->assembleMatrix("K");
}

Matrix<Real> ContactMechanicsInternodesModel::testConstructInterfaceMass(NodeGroup & contact_node_group) {
  return constructInterfaceMass(contact_node_group);
}

} // namespace akantu
