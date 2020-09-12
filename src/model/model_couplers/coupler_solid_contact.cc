/**
 * @file   coupler_solid_contact.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Thu Jan 17 2019
 * @date last modification: Thu May 22 2019
 *
 * @brief  class for coupling of solid mechanics and conatct mechanics
 * model
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
#include "coupler_solid_contact.hh"
#include "dumpable_inline_impl.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

CouplerSolidContact::CouplerSolidContact(
    Mesh & mesh, UInt dim, const ID & id, const MemoryID & memory_id,
    std::shared_ptr<DOFManager> dof_manager, const ModelType model_type)
    : Model(mesh, model_type, dof_manager, dim, id, memory_id) {

  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("CouplerSolidContact", mesh,
                                               Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("coupler_solid_contact", id, true);
  this->mesh.addDumpMeshToDumper("coupler_solid_contact", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_regular);
#endif

  this->registerDataAccessor(*this);

  solid =
      new SolidMechanicsModel(mesh, Model::spatial_dimension,
                              "solid_mechanics_model", 0, this->dof_manager);
  contact = new ContactMechanicsModel(mesh, Model::spatial_dimension,
                                      "contact_mechanics_model", 0,
                                      this->dof_manager);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
CouplerSolidContact::~CouplerSolidContact() {}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initFullImpl(const ModelOptions & options) {

  Model::initFullImpl(options);

  this->initBC(*this, *displacement, *displacement_increment, *external_force);

  solid->initFull(  _analysis_method  = this->method); 
  contact->initFull(_analysis_method  = this->method);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initModel() {

  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
FEEngine & CouplerSolidContact::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::initSolver(TimeStepSolverType time_step_solver_type,
				     NonLinearSolverType non_linear_solver_type) {
  auto & solid_model_solver =
    aka::as_type<ModelSolver>(*solid);
  solid_model_solver.initSolver(time_step_solver_type,  non_linear_solver_type);

  auto & contact_model_solver =
    aka::as_type<ModelSolver>(*contact);
  contact_model_solver.initSolver(time_step_solver_type,  non_linear_solver_type);
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
CouplerSolidContact::getDefaultSolverID(const AnalysisMethod & method) {

  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
  }
  case _explicit_consistent_mass: {
    return std::make_tuple("explicit", TimeStepSolverType::_dynamic);
  }
  case _static: {
    return std::make_tuple("static", TimeStepSolverType::_static);
  }
  case _implicit_dynamic: {
    return std::make_tuple("implicit", TimeStepSolverType::_dynamic);
  }
  default:
    return std::make_tuple("unknown", TimeStepSolverType::_not_defined);
  }  
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType CouplerSolidContact::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions CouplerSolidContact::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type =
        NonLinearSolverType::_newton_raphson_contact;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
    break;
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleResidual() {

  // computes the internal forces
  
  switch (method) {
  case _explicit_lumped_mass: {    
    auto & current_positions = contact->getContactDetector().getPositions();
    current_positions.copy(solid->getCurrentPosition());
    contact->search();
    break;
  }
  default:
    break;
  }

  this->assembleInternalForces();

  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();

  auto & contact_force = contact->getInternalForce();

  

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement", external_force, 1);
  this->getDOFManager().assembleToResidual("displacement", internal_force, 1);
  this->getDOFManager().assembleToResidual("displacement", contact_force, 1);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  //contact->assembleInternalForces();
  
  auto & internal_force = solid->getInternalForce();
  auto & external_force = solid->getExternalForce();
  auto & contact_force = contact->getInternalForce();
  
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
void CouplerSolidContact::predictor() {
  
  auto & solid_model_solver =
    aka::as_type<ModelSolver>(*solid);
  solid_model_solver.predictor(); 
      
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::corrector() {

  auto & solid_model_solver =
    aka::as_type<ModelSolver>(*solid);
  solid_model_solver.corrector();
  
  switch (method) {
  case _static:
  case _implicit_dynamic:  {
    auto & current_positions = contact->getContactDetector().getPositions();
    current_positions.copy(solid->getCurrentPosition());
    contact->search();
    break;
  }
  default:
    break;
  }

}

/* -------------------------------------------------------------------------- */
MatrixType CouplerSolidContact::getMatrixType(const ID & matrix_id) {

  if (matrix_id == "K")
    return _symmetric;
  if (matrix_id == "M") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    solid->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::beforeSolveStep() {
  auto & solid_solver_callback =
    aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.beforeSolveStep();
  
  auto & contact_solver_callback =
    aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::afterSolveStep(bool converged) {
  auto & solid_solver_callback =
    aka::as_type<SolverCallback>(*solid);
  solid_solver_callback.afterSolveStep(converged);
  
  auto & contact_solver_callback =
    aka::as_type<SolverCallback>(*contact);
  contact_solver_callback.afterSolveStep(converged);
}


  
/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  solid->assembleInternalForces();
  contact->assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  solid->assembleStiffnessMatrix(true);
 
  switch (method) {
  case _static:
  case _implicit_dynamic: {
    contact->assembleStiffnessMatrix();    
    break;
  }
  default:
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMassLumped() { solid->assembleMassLumped(); }

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMass() { solid->assembleMass(); }

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMassLumped(GhostType ghost_type) {
  solid->assembleMassLumped(ghost_type);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::assembleMass(GhostType ghost_type) {
  solid->assembleMass(ghost_type);
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> CouplerSolidContact::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    const ElementKind & kind) {

  std::shared_ptr<dumpers::Field> field;
  field = solid->createElementalField(field_name, group_name, padding_flag,
                                     spatial_dimension, kind);
 
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
CouplerSolidContact::createNodalFieldReal(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  std::shared_ptr<dumpers::Field> field;
  if (field_name == "contact_force" or field_name == "normals" or
      field_name == "normal_force" or field_name == "tangential_force" or
      field_name == "contact_state" or
      field_name == "gaps" or field_name == "previous_gaps" or
      field_name == "areas" or field_name == "tangents") {
    field =  contact->createNodalFieldReal(field_name, group_name, padding_flag);
  } else {
    field = solid->createNodalFieldReal(field_name, group_name, padding_flag);
  }

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
CouplerSolidContact::createNodalFieldBool(const std::string & field_name,
                                          const std::string & group_name,
                                          bool padding_flag) {

  
  std::shared_ptr<dumpers::Field> field;
  field = solid->createNodalFieldBool(field_name, group_name, padding_flag);
  return field;
}

#else

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
CouplerSolidContact::createElementalField(const std::string &,
                                          const std::string &, bool,
                                          const UInt &, const ElementKind &) {
  return nullptr;
}

/* ----------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field>
CouplerSolidContact::createNodalFieldReal(const std::string &,
                                          const std::string &, bool) {
  return nullptr;
}

/*-------------------------------------------------------------------*/
std::shared_ptr<dumpers::Field>
CouplerSolidContact::createNodalFieldBool(const std::string &,
                                          const std::string &, bool) {
  return nullptr;
}

#endif

/* --------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name) {
  solid->onDump();
  mesh.dump(dumper_name);
}

/* --------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name, UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, step);
}

/* -------------------------------------------------------------------------
 */
void CouplerSolidContact::dump(const std::string & dumper_name, Real time,
                               UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump() {
  solid->onDump();
  mesh.dump();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump(UInt step) {
  solid->onDump();
  mesh.dump(step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidContact::dump(Real time, UInt step) {
  solid->onDump();
  mesh.dump(time, step);
}

} // namespace akantu
