/**
 * @file   poisson_model.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Srinivasa Babu Ramisetti <srinivasa.ramisetti@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Rui Wang <rui.wang@epfl.ch>
 *
 * @date creation: Sun May 01 2011
 * @date last modification: Fri Apr 09 2021
 *
 * @brief  Implementation of PoissonModel class
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
#include "poisson_model.hh"
#include "dumpable_inline_impl.hh"
#include "element_synchronizer.hh"
#include "fe_engine_template.hh"
#include "generalized_trapezoidal.hh"
#include "group_manager_inline_impl.hh"
#include "integrator_gauss.hh"
#include "mesh.hh"
#include "parser.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include "dumper_element_partition.hh"
#include "dumper_elemental_field.hh"
#include "dumper_internal_material_field.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */



namespace akantu {
  
  class ComputeEffectiveCapacityFunctor {
  public:
    ComputeEffectiveCapacityFunctor(const PoissonModel & model) : model(model){};
    
    void operator()(Matrix<Real> & rho, const Element & element) {
      const Array<UInt> & law_indexes =
        model.getConstitutiveLawByElement(element.type, element.ghost_type);
      Real law_rho = 
        model.getConstitutiveLaw(law_indexes(element.element)).getEffectiveCapacity();
      rho.set(law_rho);
    }
    
  private:
    const PoissonModel & model;
  };
  

/* -------------------------------------------------------------------------- */
PoissonModel::PoissonModel(Mesh & mesh, UInt dim, const ID & id,
			   std::shared_ptr<DOFManager> dof_manager)
  : Model(mesh, ModelType::_poisson_model, std::move(dof_manager), dim, id),
      constitutive_law_index("constitutive law index", id),
      constitutive_law_local_numbering("constitutive law local numbering", id) {
  
  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<FEEngineType>("PoissonModelFEEngine", mesh,
					     Model::spatial_dimension);
   
  this->mesh.registerDumper<DumperParaview>("poisson_model", id, true);
  this->mesh.addDumpMesh(mesh, Model::spatial_dimension, _not_ghost,
			 _ek_regular);
  
  constitutive_law_selector =
      std::make_shared<DefaultConstitutiveLawSelector>(constitutive_law_index);

  this->registerDataAccessor(*this);

  if (this->mesh.isDistributed()) {
    auto & synchronizer = this->mesh.getElementSynchronizer();
    this->registerSynchronizer(synchronizer, SynchronizationTag::_constitutive_law_id);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_pm_dof);
    this->registerSynchronizer(synchronizer, SynchronizationTag::_pm_gradient_dof);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
PoissonModel::~PoissonModel() = default;


/* -------------------------------------------------------------------------- */
void PoissonModel::initFullImpl(const ModelOptions & options) {
  Model::initFullImpl(options);

  constitutive_law_index.initialize(mesh, _element_kind = _ek_not_defined,
				    _default_value = UInt(-1), _with_nb_element = true);
  constitutive_law_local_numbering.initialize(mesh, _element_kind = _ek_not_defined,
					      _with_nb_element = true);

  Model::initFullImpl(options);

  // initialize the materials
  if (not this->parser.getLastParsedFile().empty()) {
    this->instantiateConstitutiveLaws();
    this->initConstitutiveLaws();
  }

  this->initBC(*this, *dof, *increment, *external_dof_rate);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType PoissonModel::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions PoissonModel::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["dof"] =
        IntegrationSchemeType::_forward_euler;
    options.solution_type["dof"] = IntegrationScheme::_temperature_rate;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
    options.integration_scheme_type["dof"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["dof"] = IntegrationScheme::_not_defined;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    if (this->method == _explicit_consistent_mass) {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["dof"] =
          IntegrationSchemeType::_forward_euler;
      options.solution_type["dof"] =
          IntegrationScheme::_temperature_rate;
    } else {
      options.non_linear_solver_type = NonLinearSolverType::_newton_raphson;
      options.integration_scheme_type["dof"] =
          IntegrationSchemeType::_backward_euler;
      options.solution_type["dof"] = IntegrationScheme::_temperature;
    }
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
  }

  return options;
}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
PoissonModel::getDefaultSolverID(const AnalysisMethod & method) {
  switch (method) {
  case _explicit_lumped_mass: {
    return std::make_tuple("explicit_lumped",
                           TimeStepSolverType::_dynamic_lumped);
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
void PoissonModel::initSolver(TimeStepSolverType time_step_solver_type,
                                   NonLinearSolverType /*unused*/) {
  auto & dof_manager = this->getDOFManager();

  this->allocNodalField(this->dof, 1, "dof");
  this->allocNodalField(this->external_dof_rate, 1, "external_dof_rate");
  this->allocNodalField(this->internal_dof_rate, 1, "internal_dof_rate");
  this->allocNodalField(this->blocked_dofs, 1, "blocked_dofs");

  if (!dof_manager.hasDOFs("dof")) {
    dof_manager.registerDOFs("dof", *this->dof, _dst_nodal);
    dof_manager.registerBlockedDOFs("dof", *this->blocked_dofs);
  }

  if (time_step_solver_type == TimeStepSolverType::_dynamic ||
      time_step_solver_type == TimeStepSolverType::_dynamic_lumped) {
    this->allocNodalField(this->dof_rate, 1, "dof_rate");

    if (!dof_manager.hasDOFsDerivatives("dof", 1)) {
      dof_manager.registerDOFsDerivative("dof", 1,
                                         *this->dof_rate);
    }
  }
}

  
/* -------------------------------------------------------------------------- */
void PoissonModel::initModel() {
  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);

}

  
/* -------------------------------------------------------------------------- */
void PoissonModel::assembleResidual() {
  AKANTU_DEBUG_IN();

  this->assembleInternalDofRate();

  this->getDOFManager().assembleToResidual("dof",
                                           *this->external_dof_rate, 1);
  this->getDOFManager().assembleToResidual("dof",
                                           *this->internal_dof_rate, 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
MatrixType PoissonModel::getMatrixType(const ID & matrix_id) const {
  if (matrix_id == "K" or matrix_id == "M") {
    return _symmetric;
  }

  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleMatrix(const ID & matrix_id) {
  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacity();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M" and need_to_reassemble_capacity) {
    this->assembleCapacityLumped();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::beforeSolveStep() {
  for (auto & law : constitutive_laws) {
    law->beforeSolveStep();
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::afterSolveStep(bool converged) {
  for (auto & law : constitutive_laws) {
    law->afterSolveStep(converged);
  }
}

/* -------------------------------------------------------------------------- */
void PoissonModel::predictor() { ++dof_release; }
 
/* -------------------------------------------------------------------------- */
void PoissonModel::corrector() { ++dof_release; }

/* -------------------------------------------------------------------------- */
void PoissonModel::assembleInternalDofRate() {
  AKANTU_DEBUG_IN();

  this->internal_dof_rate->zero();
  
  this->synchronize(SynchronizationTag::_pm_dof);

  // compute the fluxes of local elements
  AKANTU_DEBUG_INFO("Compute local fluxes");
  for (auto & law : constitutive_laws) {
    law->computeAllFluxes(_not_ghost);
  }

 
  // assemble the rate due to local fluxes
  AKANTU_DEBUG_INFO("Assemble residual for local elements");
  for (auto & law : constitutive_laws) {
    law->assembleInternalDofRate(_not_ghost);
  }

  // compute the fluxes of ghost elements
  AKANTU_DEBUG_INFO("Compute local fluxes");
  for (auto & law : constitutive_laws) {
    law->computeAllFluxes(_ghost);
  }
  
  // assemble the fluxes due to ghost elements
  AKANTU_DEBUG_INFO("Assemble residual for ghost elements");
  for (auto & law : constitutive_laws) {
    law->assembleInternalDofRate(_ghost);
  }

 
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleStiffnessMatrix(bool need_to_reassemble) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix.");

  // Check if constitutive laws need to recompute the matrix
  for (auto & law : constitutive_laws) {
    need_to_reassemble |= law->hasMatrixChanged("K");
  }

  if (need_to_reassemble) {
    this->getDOFManager().getMatrix("K").zero();

    // call compute stiffness matrix on each local elements
    for (auto & law : constitutive_laws) {
      law->assembleStiffnessMatrix(_not_ghost);
    }
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacityLumped() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_capacity_lumped) {
    return;
  }

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  assembleCapacityLumped(_not_ghost);
  assembleCapacityLumped(_ghost);
    
  need_to_reassemble_capacity_lumped = false;

  AKANTU_DEBUG_OUT();
}

  
/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacityLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<FEEngineType>();
  ComputeEffectiveCapacityFunctor compute_rho(*this);

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldLumped(compute_rho, "M", "dof",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
void PoissonModel::assembleCapacity() {
  AKANTU_DEBUG_IN();
  auto ghost_type = _not_ghost;

  this->getDOFManager().zeroMatrix("M");

  auto & fem = getFEEngineClass<FEEngineType>();

  ComputeEffectiveCapacityFunctor rho_functor(*this);

  for (auto && type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(rho_functor, "M", "dof",
                            this->getDOFManager(), type, ghost_type);
  }

  need_to_reassemble_capacity = false;

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
FEEngine & PoissonModel::getFEEngineBoundary(const ID & name) {
  return aka::as_type<FEEngine>(getFEEngineClassBoundary<FEEngineType>(name));
}


/* -------------------------------------------------------------------------- */
/* Information                                                                */
/* -------------------------------------------------------------------------- */
Real PoissonModel::getStableTimeStep() {
  AKANTU_DEBUG_IN();

  Real min_dt = getStableTimeStep(_not_ghost);

  /// reduction min over all processors
  mesh.getCommunicator().allReduce(min_dt, SynchronizerOperation::_min);

  AKANTU_DEBUG_OUT();
  return min_dt;
}

/* -------------------------------------------------------------------------- */
Real PoissonModel::getStableTimeStep(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real min_dt = std::numeric_limits<Real>::max();

  Element elem;
  elem.ghost_type = ghost_type;

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    elem.type = type;
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
    UInt nb_element = mesh.getNbElement(type);

    auto law_indexes = constitutive_law_index(type, ghost_type).begin();
    auto law_loc_num = constitutive_law_local_numbering(type, ghost_type).begin();

    Array<Real> X(0, nb_nodes_per_element * Model::spatial_dimension);
    FEEngine::extractNodalToElementField(mesh, mesh.getNodes(), X, type,
                                         _not_ghost);

    auto X_el = X.begin(Model::spatial_dimension, nb_nodes_per_element);

    for (UInt el = 0; el < nb_element;
         ++el, ++X_el, ++law_indexes, ++law_loc_num) {
      elem.element = *law_loc_num;
      Real el_h = getFEEngine().getElementInradius(*X_el, type);
      Real el_c = this->constitutive_laws[*law_indexes]->getCelerity();
      Real el_dt = el_h / el_c;

      min_dt = std::min(min_dt, el_dt);
    }
  }

  AKANTU_DEBUG_OUT();
  return min_dt;
}
  

/* -------------------------------------------------------------------------- */

void PoissonModel::setTimeStep(Real time_step, const ID & solver_id) {
  Model::setTimeStep(time_step, solver_id);

  this->mesh.getDumper("poisson_model").setTimeStep(time_step);
}


/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createNodalFieldBool(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  std::map<std::string, Array<bool> *> uint_nodal_fields;
  uint_nodal_fields["blocked_dofs"] = blocked_dofs.get();

  auto field = mesh.createNodalField(uint_nodal_fields[field_name], group_name);
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createNodalFieldReal(
    const std::string & field_name, const std::string & group_name,
    __attribute__((unused)) bool padding_flag) {

  if (field_name == "capacity_lumped") {
    AKANTU_EXCEPTION(
        "Capacity lumped is a nodal field now stored in the DOF manager."
        "Therefore it cannot be used by a dumper anymore");
  }

  std::map<std::string, Array<Real> *> real_nodal_fields;
  real_nodal_fields["dof"] = dof.get();
  real_nodal_fields["dof_rate"] = dof_rate.get();
  real_nodal_fields["external_dof_rate"] = external_dof_rate.get();
  real_nodal_fields["internal_dof_rate"] = internal_dof_rate.get();
  real_nodal_fields["increment"] = increment.get();

  std::shared_ptr<dumpers::Field> field =
      mesh.createNodalField(real_nodal_fields[field_name], group_name);

  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumpers::Field> PoissonModel::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool /*padding_flag*/, UInt /*spatial_dimension*/,
    ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;

  if (field_name == "partitions") {
    field = mesh.createElementalField<UInt, dumpers::ElementPartitionField>(
        mesh.getConnectivities(), group_name, this->spatial_dimension,
        element_kind);
  }
  
  return field;
}

/* -------------------------------------------------------------------------- */
inline UInt PoissonModel::getNbData(const Array<UInt> & indexes,
				    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  
  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::packData(CommunicationBuffer & buffer,
                                        const Array<UInt> & indexes,
                                        const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::unpackData(CommunicationBuffer & buffer,
                                          const Array<UInt> & indexes,
                                          const SynchronizationTag & tag) {
  AKANTU_DEBUG_IN();


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
inline UInt PoissonModel::getNbData(const Array<Element> & elements,
				    const SynchronizationTag & tag) const {
  AKANTU_DEBUG_IN();

  UInt size = 0;
  UInt nb_nodes_per_element = 0;

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::packData(CommunicationBuffer & buffer,
				   const Array<Element> & elements,
				   const SynchronizationTag & tag) const {
  
}

/* -------------------------------------------------------------------------- */
inline void PoissonModel::unpackData(CommunicationBuffer & buffer,
				     const Array<Element> & elements,
				     const SynchronizationTag & tag) {
  
}

/* -------------------------------------------------------------------------- */
} // namespace akantu
