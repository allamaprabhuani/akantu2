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
      dof_manager);
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

  // init solid mechanics model
  solid->initFull(options);

  detector->findContactNodes();

  Model::initFullImpl(options);
}
/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleResidual() {
  // get positions of master and slave nodes
  auto && positions = mesh.getNodes();
  auto && positions_view = make_view(positions, spatial_dimension).begin();

  auto & master_node_group = detector->getMasterNodeGroup();
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
  Array<Real> constraints(master_node_group.size(), spatial_dimension);

  for (int dim = 0; dim < spatial_dimension; dim++) {
    auto tmp = R_master_slave * positions_slave(dim);
    for (UInt i = 0; i < master_node_group.size(); i++) {
      constraints(i, dim) = tmp[i] - positions_master(i, dim);
    }
  }
  dof_manager->assembleToResidual("lambdas", constraints, 1);
}

MatrixType
ContactMechanicsInternodesModel::getMatrixType(const ID & matrix_id) const {
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    solid->assembleStiffnessMatrix();

    auto & master_node_group = detector->getMasterNodeGroup();
    auto & slave_node_group = detector->getSlaveNodeGroup();

    auto nb_master_nodes = master_node_group.getNodes().size();
    auto nb_slave_nodes = slave_node_group.getNodes().size();

    auto nb_nodes = mesh.getNbNodes();
    auto nb_dofs = spatial_dimension * nb_nodes;
    auto nb_constraint_dofs = spatial_dimension * nb_master_nodes;

    // R matrices, need for calculations of blocks
    auto && R_master_slave = detector->constructInterpolationMatrix(
        master_node_group, slave_node_group, detector->getSlaveRadiuses());
    auto && R_slave_master = detector->constructInterpolationMatrix(
        slave_node_group, master_node_group, detector->getMasterRadiuses());

    // assemble B
    TermsToAssemble termsB("displacement", "lambdas");
    for (UInt i = 0; i < nb_dofs; ++i) {
      for (UInt j = 0; j < nb_constraint_dofs; ++j) {
        termsB(i, j) = 1;
      }
    }

    // assemble B_tilde
    TermsToAssemble termsB_tilde("lambdas", "displacement");
    for (UInt i = 0; i < nb_constraint_dofs; ++i) {
      for (UInt j = 0; j < nb_dofs; ++j) {
        termsB_tilde(i, j) = 1.;
      }
    }

    // assemble C (zeros)
    TermsToAssemble termsC("lambdas", "lambdas");
    for (UInt i = 0; i < nb_constraint_dofs; ++i) {
      for (UInt j = 0; j < nb_constraint_dofs; ++j) {
        termsB_tilde(i, j) = 1.;
      }
    }

    dof_manager->assemblePreassembledMatrix("K", termsB);
    dof_manager->assemblePreassembledMatrix("K", termsB_tilde);
    dof_manager->assemblePreassembledMatrix("K", termsC);

    auto & A = dof_manager->getMatrix("K");
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
  // not sure how to initalize initSolver
  solid->solveStep(callback, solver_id);

  ModelSolver::solveStep(callback, solver_id);
}

void ContactMechanicsInternodesModel::solveStep(const ID & solver_id) {

  ModelSolver::solveStep(solver_id);
}

/* -------------------------------------------------------------------------- */
void ContactMechanicsInternodesModel::beforeSolveStep() {}

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
    TimeStepSolverType /*unused*/, NonLinearSolverType /*unused*/) {

  auto nb_master_nodes = detector->getMasterNodeGroup().getNodes().size();

  // allocate lambdas
  lambdas = std::make_unique<Array<Real>>(nb_master_nodes, spatial_dimension,
                                          id + ":lambdas");

  // allocate blocked dofs for lambdas
  blocked_dofs = std::make_unique<Array<bool>>(
      nb_master_nodes, spatial_dimension, id + ":blocked_dofs");

  dof_manager->registerDOFs("lambdas", *lambdas, _dst_generic);
  dof_manager->registerBlockedDOFs("lambdas", *blocked_dofs);
}

} // namespace akantu
