/**
 * @file   coupler_solid_contact_tmpl.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jan 21 2019
 * @date last modification: Wed Jun 23 2021
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "coupler_solid_contact.hh"
#include "dumpable_inline_impl.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::AbstractCouplerSolidContactTemplate(Mesh & mesh, const ModelType & type,
                                          UInt dim, const ID & id,
                                          std::shared_ptr<DOFManager> dof_manager)
    : Model(mesh, type, dof_manager, dim, id) {
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::initSolver(
    TimeStepSolverType time_step_solver_type,
    NonLinearSolverType non_linear_solver_type) {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.initSolver(time_step_solver_type, non_linear_solver_type);

  auto & contact_model_solver = aka::as_type<ModelSolver>(*contact);
  contact_model_solver.initSolver(time_step_solver_type,
                                  non_linear_solver_type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
std::tuple<ID, TimeStepSolverType>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::getDefaultSolverID(
    const AnalysisMethod & method) {
  return solid->getDefaultSolverID(method);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
TimeStepSolverType
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::getDefaultSolverType()
    const {
  return solid->getDefaultSolverType();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
ModelSolverOptions
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  return solid->getDefaultSolverOptions(type);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::predictor() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.predictor();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::corrector() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(*solid);
  solid_model_solver.corrector();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::assembleLumpedMatrix(
    const ID & matrix_id) {
  if (matrix_id == "M") {
    solid->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::beforeSolveStep() {
  auto & solid_solver_callback = aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.beforeSolveStep();

  auto & contact_solver_callback = aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::afterSolveStep(bool converged) {
  auto & solid_solver_callback = aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.afterSolveStep(converged);

  auto & contact_solver_callback = aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.afterSolveStep(converged);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
std::shared_ptr<dumpers::Field>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, UInt spatial_dimension, ElementKind kind) {

  std::shared_ptr<dumpers::Field> field;
  field = contact->createElementalField(field_name, group_name, padding_flag,
                                        spatial_dimension, kind);
  if (not field) {
    field = solid->createElementalField(field_name, group_name, padding_flag,
                                        spatial_dimension, kind);
  }
  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
std::shared_ptr<dumpers::Field>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {

  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldReal(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldReal(field_name, group_name, padding_flag);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
std::shared_ptr<dumpers::Field>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>::
    createNodalFieldUInt(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {
  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldUInt(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldUInt(field_name, group_name, padding_flag);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
std::shared_ptr<dumpers::Field>
AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag) {
  std::shared_ptr<dumpers::Field> field;
  field = contact->createNodalFieldBool(field_name, group_name, padding_flag);
  if (not field) {
    field = solid->createNodalFieldBool(field_name, group_name, padding_flag);
  }
  return field;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump(const std::string & dumper_name) {
  solid->onDump();
  Model::dump(dumper_name);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump(const std::string & dumper_name, UInt step) {
  solid->onDump();
  Model::dump(dumper_name, step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump(const std::string & dumper_name, Real time, UInt step) {
  solid->onDump();
  Model::dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump() {
  solid->onDump();
  Model::dump();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump(UInt step) {
  solid->onDump();
  Model::dump(step);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType, class ContactMechanicsModelType>
void AbstractCouplerSolidContactTemplate<SolidMechanicsModelType, ContactMechanicsModelType>
    ::dump(Real time, UInt step) {
  solid->onDump();
  Model::dump(time, step);
}

/* -------------------------------------------------------------------------- */
/* Implementations specific to ContactMechanicsModel                          */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
MatrixType CouplerSolidContactTemplate<SolidMechanicsModelType>::getMatrixType(
    const ID & matrix_id) const {
  if (matrix_id == "K") {
    return _symmetric;
  }
  if (matrix_id == "M") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleMatrix(
        const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    this->solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::corrector() {
  AbstractCouplerSolidContactTemplate<SolidMechanicsModelType,
                                      ContactMechanicsModel>::corrector();

  switch (this->method) {
  case _static:
  case _implicit_dynamic: {
    auto & current_positions = this->contact->getContactDetector().getPositions();
    current_positions.copy(this->solid->getCurrentPosition());
    this->contact->search();
    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleResidual() {
  // computes the internal forces
  switch (this->method) {
  case _explicit_lumped_mass: {
    auto & current_positions = this->contact->getContactDetector().getPositions();
    current_positions.copy(this->solid->getCurrentPosition());
    this->contact->search();
    break;
  }
  default:
    break;
  }

  this->assembleInternalForces();

  auto & internal_force = this->solid->getInternalForce();
  auto & external_force = this->solid->getExternalForce();

  auto & contact_force = this->contact->getInternalForce();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement", external_force, 1);
  this->getDOFManager().assembleToResidual("displacement", internal_force, 1);
  this->getDOFManager().assembleToResidual("displacement", contact_force, 1);
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleResidual(
        const ID & residual_part) {
  AKANTU_DEBUG_IN();

  // contact->assembleInternalForces();

  auto & internal_force = this->solid->getInternalForce();
  auto & external_force = this->solid->getExternalForce();
  auto & contact_force = this->contact->getInternalForce();

  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement", external_force, 1);
    this->getDOFManager().assembleToResidual("displacement", contact_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement", internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  this->solid->assembleInternalForces();
  this->contact->assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class SolidMechanicsModelType>
void CouplerSolidContactTemplate<SolidMechanicsModelType>::
    assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  this->solid->assembleStiffnessMatrix(true);

  switch (this->method) {
  case _static:
  case _implicit_dynamic: {
    this->contact->assembleStiffnessMatrix();
    break;
  }
  default:
    break;
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
